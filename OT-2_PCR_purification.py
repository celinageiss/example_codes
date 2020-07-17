# =============================================================================
# PCR purification protocol
# =============================================================================

# Define the parameters

def get_values(*names):
    import json
    _all_values = json.loads("""{
        "pipette_type": "p300_multi",
        "pipette_mount": "right",
        "column_number": 12,
        "PCR_volume": 20,
        "bead_ratio": 1.8,
        "elution_buffer_volume": 60,
        "incubation_time": 5,
        "settling_time": 60,
        "drying_time": 15
        }""")
    return [_all_values[n] for n in names]

# =============================================================================
# Now load the protocol on the robot. Please do NOT change the script below!
# =============================================================================


import math

metadata = {
    'protocolName': 'PCR purification using magnetic beads (Omega)',
    'author': 'Celina Geiss <celina.geiss@dkfz-heidelberg.de>',
    'apiLevel': '2.1'
}


def run(protocol):
    [pipette_type, pipette_mount, column_number, PCR_volume, bead_ratio,
     elution_buffer_volume, incubation_time, settling_time,
     drying_time] = get_values(  # noqa: F821
        "pipette_type", "pipette_mount", "column_number", "PCR_volume",
        "bead_ratio", "elution_buffer_volume", "incubation_time",
        "settling_time", "drying_time")

    mag_deck = protocol.load_module('magdeck', 1)
    mag_plate = mag_deck.load_labware(
        'biorad_96_wellplate_200ul_pcr', label='PCR plate')
    reagent_container = protocol.load_labware(
        'usascientific_12_reservoir_22ml', 2, label='reagent reservoir')
    output_plate = protocol.load_labware(
        'biorad_96_wellplate_200ul_pcr', 3, label='elution plate')
    liquid_waste = protocol.load_labware('agilent_1_reservoir_290ml', 9, label='liquid waste')
    liquid_waste_top = liquid_waste.wells()[0].top()

    # Define tips and pipette
    sample_number = column_number*8
    total_tips = sample_number*8
    tiprack_num = math.ceil(total_tips/96)
    tipslots = [4, 5, 7, 8, 10, 11][:tiprack_num]

    tipracks = [protocol.load_labware('opentrons_96_tiprack_300ul', slot,
                                      label='Neptune 200 uL tiprack') for slot in tipslots]

    pipette = protocol.load_instrument(
        pipette_type, pipette_mount, tip_racks=tipracks)

    col_num = math.ceil(sample_number/8)
    samples = [col for col in mag_plate.rows()[0][:col_num]]
    samples_top = [well.top() for well in mag_plate.rows()[0][:col_num]]
    output = [col for col in output_plate.rows()[0][:col_num]]

    # Define reagents and liquid waste
    beads = reagent_container.wells()[0]            # A1

    ethanol_1 = reagent_container.wells()[2]        # A3
    ethanol_2 = reagent_container.wells()[3]        # A4
    ethanol_3 = reagent_container.wells()[4]        # A5
    ethanol_4 = reagent_container.wells()[5]        # A6
    ethanol_wells = [ethanol_1, ethanol_2, ethanol_3, ethanol_4]

    elution_buffer = reagent_container.wells()[7]   # A8

    # Disengage MagDeck
    mag_deck.disengage()

    # Define bead and mix volume to resuspend beads
    bead_volume = PCR_volume*bead_ratio
    total_vol = bead_volume + PCR_volume + 15
    mix_vol_target = total_vol/2

    # Mix beads and PCR samples
    protocol.comment("Resuspending beads.")
    pipette.flow_rate.aspirate = 180
    pipette.flow_rate.dispense = 180
    pipette.pick_up_tip()
    pipette.mix(15, 200, beads.bottom(z=3))

    protocol.comment("Mixing beads and PCR samples.")
    mix_vol = 200
    for target in samples:
        if not pipette.hw_pipette['has_tip']:
            pipette.pick_up_tip()
        pipette.flow_rate.aspirate = 180
        pipette.flow_rate.dispense = 180
        pipette.mix(3, mix_vol, beads)
        protocol.default_speed = 200  # Slow down head speed 0.5X for bead handling
        pipette.flow_rate.aspirate = 10
        pipette.flow_rate.dispense = 10
        pipette.transfer(bead_volume, beads, target, new_tip='never')
        pipette.flow_rate.aspirate = 50
        pipette.flow_rate.dispense = 50
        pipette.mix(25, mix_vol_target, target)  # originally 40
        pipette.blow_out()
        protocol.default_speed = 400
        pipette.drop_tip()

    # Incubate beads and PCR product at RT for 5 minutes
    protocol.comment("Incubating the beads and PCR products at room \
temperature for %d minutes. Protocol will resume automatically." % incubation_time)
    protocol.delay(minutes=incubation_time)

    # Engage MagDeck and Magnetize
    mag_deck.engage(height=19)
    protocol.delay(seconds=settling_time)

    # Remove supernatant from magnetic beads
    pipette.flow_rate.aspirate = 25
    pipette.flow_rate.dispense = 120
    for target in samples:
        pipette.transfer(
            total_vol + 10, target, liquid_waste_top, blow_out=True)

    # Wash beads twice with 70% ethanol
    air_vol = 15

    for _ in range(2):
        pipette.pick_up_tip()
        for i in range(0, column_number):
            ethanol = ethanol_wells[(i//3)]  # take new EtOH well every third column
            pipette.transfer(
                160, ethanol, samples_top[i], air_gap=air_vol, new_tip='never')
        protocol.delay(minutes=1)
        for target in samples:
            if not pipette.hw_pipette['has_tip']:
                pipette.pick_up_tip()
            pipette.transfer(180, target.bottom(z=0.5), liquid_waste_top,
                             air_gap=air_vol, blow_out=True, new_tip='never')
        pipette.drop_tip()

    # Dry at RT
    msg = "Drying the beads for %d minutes. Protocol will resume automatically." % drying_time
    protocol.delay(minutes=drying_time, msg=msg)

    # Disengage MagDeck
    mag_deck.disengage()

    # Mix beads with elution buffer
    mix_vol = elution_buffer_volume/2
    for target in samples:
        pipette.transfer(elution_buffer_volume, elution_buffer, target,
                         mix_after=(45, mix_vol))

    # Incubate at RT for 3 minutes
    protocol.comment("Incubating at room temperature for 3 minutes. \
Protocol will resume automatically.")
    protocol.delay(minutes=3)

    # Engage MagDeck for settling_time and remain engaged for DNA elution
    mag_deck.engage(height=22)
    protocol.comment("Delaying for %d seconds for \
beads to settle." % settling_time)
    protocol.delay(seconds=settling_time)

    # Transfer clean PCR product to a new well
    pipette.flow_rate.aspirate = 25
    pipette.flow_rate.dispense = 120
    for target, dest in zip(samples, output):
        pipette.pick_up_tip()
        pipette.transfer(elution_buffer_volume, target.bottom(z=1), dest.top(), new_tip='never')
        pipette.blow_out(dest.top(z=-2))
        pipette.drop_tip()

    # Disengage MagDeck
    mag_deck.disengage()

    # End of protocol, home robot
    protocol.home()
    protocol.comment("DNA purification protocol finished. Please take out your purified samples (slot 2).")

