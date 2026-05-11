from rt_qpcr_guide.rt_calculator import calculate_rt_reaction, plan_reactions


def test_calculate_rt_reaction_uses_equal_rna_input():
    reaction = calculate_rt_reaction(
        sample_id="S1",
        concentration_ng_ul=250,
        target_rna_ng=1000,
    )

    assert reaction.rna_volume_ul == 4.0
    assert reaction.gdna_wiper_mix_ul == 4.0
    assert reaction.water_to_16_ul == 8.0
    assert reaction.rt_supermix_ul == 4.0
    assert reaction.total_volume_ul == 20.0
    assert reaction.status == "ok"


def test_calculate_rt_reaction_flags_low_concentration():
    reaction = calculate_rt_reaction(
        sample_id="Low",
        concentration_ng_ul=80,
        target_rna_ng=1000,
    )

    assert reaction.rna_volume_ul == 12.5
    assert reaction.water_to_16_ul == 0
    assert reaction.status == "rna_volume_too_high"


def test_plan_reactions_parses_concentration_rows():
    rows = [
        {"sample_id": "S1", "concentration_ng_ul": "250"},
        {"sample_id": "S2", "concentration_ng_ul": "200"},
    ]

    plan = plan_reactions(rows, target_rna_ng=1000)

    assert [r.sample_id for r in plan] == ["S1", "S2"]
    assert [r.rna_volume_ul for r in plan] == [4.0, 5.0]

