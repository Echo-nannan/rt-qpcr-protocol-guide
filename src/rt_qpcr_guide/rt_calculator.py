"""Reverse transcription reaction planning."""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class RTReaction:
    sample_id: str
    concentration_ng_ul: float
    target_rna_ng: float
    rna_volume_ul: float
    gdna_wiper_mix_ul: float
    water_to_16_ul: float
    rt_supermix_ul: float
    total_volume_ul: float
    status: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def calculate_rt_reaction(
    sample_id: str,
    concentration_ng_ul: float,
    target_rna_ng: float = 1000,
    gdna_wiper_mix_ul: float = 4,
    first_stage_volume_ul: float = 16,
    rt_supermix_ul: float = 4,
) -> RTReaction:
    if concentration_ng_ul <= 0:
        raise ValueError("concentration_ng_ul must be positive")
    if target_rna_ng <= 0:
        raise ValueError("target_rna_ng must be positive")

    rna_volume = target_rna_ng / concentration_ng_ul
    water = first_stage_volume_ul - gdna_wiper_mix_ul - rna_volume
    status = "ok" if water >= 0 else "rna_volume_too_high"

    return RTReaction(
        sample_id=sample_id,
        concentration_ng_ul=round(concentration_ng_ul, 4),
        target_rna_ng=round(target_rna_ng, 4),
        rna_volume_ul=round(rna_volume, 4),
        gdna_wiper_mix_ul=round(gdna_wiper_mix_ul, 4),
        water_to_16_ul=round(max(water, 0), 4),
        rt_supermix_ul=round(rt_supermix_ul, 4),
        total_volume_ul=round(first_stage_volume_ul + rt_supermix_ul, 4),
        status=status,
    )


def plan_reactions(rows: list[dict[str, str]], target_rna_ng: float = 1000) -> list[RTReaction]:
    plan: list[RTReaction] = []
    for row in rows:
        sample_id = row.get("sample_id", "").strip()
        concentration_text = row.get("concentration_ng_ul", "").strip()
        if not sample_id or not concentration_text:
            continue
        plan.append(
            calculate_rt_reaction(
                sample_id=sample_id,
                concentration_ng_ul=float(concentration_text),
                target_rna_ng=target_rna_ng,
            )
        )
    return plan

