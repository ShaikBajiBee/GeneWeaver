"""
Sliding-window GC Skew analysis and replication origin (oriC) prediction module.
Computes local nucleotide skew (G - C)/(G + C), cumulative skew curves, and generates SVG profiles.
"""

from typing import Dict, List, Any
from analyzer import clean_sequence, SequenceError


def calculate_sliding_gc_skew(
    sequence_raw: str,
    window_size: int = 20,
    step_size: int = 5,
) -> Dict[str, Any]:
    """
    Calculate sliding-window GC skew across a DNA sequence.
    Locates the minimum of cumulative GC skew (predicted oriC / replication origin)
    and the maximum (predicted terminus ter).
    """
    seq = clean_sequence(sequence_raw)
    seq_len = len(seq)

    if seq_len < window_size:
        window_size = max(5, seq_len)
        step_size = max(1, window_size // 4)

    positions: List[int] = []
    window_skews: List[float] = []
    cumulative_skew: List[float] = []

    running_cum = 0.0

    for i in range(0, seq_len - window_size + 1, step_size):
        sub = seq[i : i + window_size]
        g = sub.count("G")
        c = sub.count("C")

        skew = (g - c) / (g + c) if (g + c) > 0 else 0.0
        skew_round = round(skew, 4)

        running_cum += skew
        midpoint = i + (window_size // 2) + 1

        positions.append(midpoint)
        window_skews.append(skew_round)
        cumulative_skew.append(round(running_cum, 4))

    # If sequence was too short for loop, take single window
    if not window_skews:
        g = seq.count("G")
        c = seq.count("C")
        sk = round((g - c) / (g + c), 4) if (g + c) > 0 else 0.0
        positions = [seq_len // 2 + 1]
        window_skews = [sk]
        cumulative_skew = [sk]

    min_cum = min(cumulative_skew)
    max_cum = max(cumulative_skew)

    ori_idx = cumulative_skew.index(min_cum)
    ter_idx = cumulative_skew.index(max_cum)

    predicted_ori = positions[ori_idx]
    predicted_ter = positions[ter_idx]

    return {
        "sequence_length": seq_len,
        "window_size": window_size,
        "step_size": step_size,
        "total_windows": len(window_skews),
        "positions": positions,
        "window_skews": window_skews,
        "cumulative_skews": cumulative_skew,
        "predicted_origin_bp": predicted_ori,
        "predicted_terminus_bp": predicted_ter,
        "min_cumulative_skew": min_cum,
        "max_cumulative_skew": max_cum,
    }


def generate_skew_svg(
    skew_data: Dict[str, Any],
    width: int = 560,
    height: int = 180,
) -> str:
    """
    Generate a lightweight inline SVG line chart visualizing cumulative GC skew.
    """
    skews = skew_data.get("cumulative_skews", [])
    if not skews:
        return ""

    min_val = min(skews)
    max_val = max(skews)
    val_range = (max_val - min_val) if max_val != min_val else 1.0

    pad_x = 40
    pad_y = 25
    plot_w = width - (2 * pad_x)
    plot_h = height - (2 * pad_y)

    points = []
    n = len(skews)
    for idx, val in enumerate(skews):
        x = pad_x + (idx / (n - 1) if n > 1 else 0.5) * plot_w
        # invert y for SVG coordinates (0 is top)
        y = pad_y + (1.0 - ((val - min_val) / val_range)) * plot_h
        points.append(f"{round(x, 1)},{round(y, 1)}")

    poly_pts = " ".join(points)
    ori_bp = skew_data.get("predicted_origin_bp", 0)

    svg = f"""<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg">
  <rect width="100%" height="100%" fill="#0f172a" rx="8" />
  <line x1="{pad_x}" y1="{height - pad_y}" x2="{width - pad_x}" y2="{height - pad_y}" stroke="#334155" stroke-width="1" />
  <polyline fill="none" stroke="#38bdf8" stroke-width="2.5" points="{poly_pts}" />
  <text x="{pad_x}" y="{pad_y - 8}" fill="#94a3b8" font-size="11" font-family="sans-serif">Cumulative GC Skew (Predicted oriC @ {ori_bp} bp)</text>
</svg>"""
    return svg
