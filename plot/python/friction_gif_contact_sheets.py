#!/usr/bin/env python3
"""Make contact sheets from friction-sweep GIF outputs."""

from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path("analysis_outputs/08jun/friction_buildup_sweep_100/plot_hd_visualizations")
OUTDIR = Path("analysis_outputs/08jun/friction_buildup_sweep_100/gif_contact_sheets")
FRICTIONS = ["friction10", "friction30", "friction100", "friction300"]
FIELDS = ["visualization", "phi", "pressure"]
FRAME_IDS = [0, 10, 20]
FRAME_LABELS = ["t=0", "t=5000", "t=10000"]


def load_font(size: int):
    for path in (
        "/System/Library/Fonts/Helvetica.ttc",
        "/System/Library/Fonts/Supplemental/Arial.ttf",
    ):
        try:
            return ImageFont.truetype(path, size)
        except OSError:
            pass
    return ImageFont.load_default()


def extract_frame(path: Path, frame_id: int, size: tuple[int, int] | None) -> Image.Image:
    gif = Image.open(path)
    gif.seek(frame_id)
    frame = gif.convert("RGB")
    return frame if size is None else frame.resize(size)


def make_sheet(field: str) -> Path:
    font = load_font(22)
    rows: list[tuple[str, Image.Image]] = []
    cell_size: tuple[int, int] | None = None

    for friction in FRICTIONS:
        gif_path = ROOT / f"{friction}_seed1001" / f"{field}_0-100_step5.gif"
        frames = []
        for frame_id in FRAME_IDS:
            frame = extract_frame(gif_path, frame_id, cell_size)
            if cell_size is None:
                cell_size = frame.size
            frames.append(frame)

        assert cell_size is not None
        row = Image.new("RGB", (cell_size[0] * len(frames), cell_size[1]), "white")
        for col, frame in enumerate(frames):
            row.paste(frame, (col * cell_size[0], 0))
        rows.append((friction, row))

    assert cell_size is not None
    label_w = 125
    header_h = 44
    sheet = Image.new(
        "RGB",
        (label_w + cell_size[0] * len(FRAME_IDS), header_h + cell_size[1] * len(rows)),
        "white",
    )
    draw = ImageDraw.Draw(sheet)
    for col, label in enumerate(FRAME_LABELS):
        draw.text((label_w + col * cell_size[0] + 10, 8), label, fill="black", font=font)

    for row_index, (friction, row) in enumerate(rows):
        y = header_h + row_index * cell_size[1]
        draw.text((8, y + cell_size[1] // 2 - 12), friction.replace("friction", "f="), fill="black", font=font)
        sheet.paste(row, (label_w, y))

    OUTDIR.mkdir(parents=True, exist_ok=True)
    outpath = OUTDIR / f"{field}_gif_contact_sheet.png"
    sheet.save(outpath)
    return outpath


def main() -> None:
    for field in FIELDS:
        print(make_sheet(field))


if __name__ == "__main__":
    main()
