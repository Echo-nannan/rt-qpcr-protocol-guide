# qPCR Result Processor Pro

A Tkinter-based qPCR data processing tool for format conversion, Delta Delta Ct analysis, and customizable exports.

## Features

### Format Conversion

- Parse raw data exported from instruments such as Roche LightCycler.
- Support 384-well plate layouts with 24 columns and 16 rows.
- Customize sample names, gene names, and group names.
- Export raw data, mean Ct values, and long-format tables for Delta Delta Ct analysis.

### Delta Delta Ct Analysis

- Run complete relative quantification with Delta Delta Ct.
- Calculate reference means by group.
- Detect available genes and sample groups automatically.
- Export detailed calculation tables and summary tables.
- Export side-by-side reference-gene and target-gene tables.

### Custom Output Names

Supported template variables:

- `{date}`: current date in `YYYYMMDD` format
- `{sample}`: user-defined sample information
- `{gene}`: gene name

Example:

```text
{date}_{sample}_Analysis and processing -> 20250124_DemoSample_Analysis and processing.csv
```

## Installation and Run

```bash
pip install -r requirements.txt
python scripts/start_gui.py
```

On Windows, `run_gui.bat` can also be used.

## Dependencies

- Python 3.8+
- pandas >= 1.5.0
- openpyxl >= 3.0.0
- numpy >= 1.21.0
- xlsxwriter >= 3.0.0
- pyyaml >= 6.0.0
- matplotlib >= 3.7.0
- scipy >= 1.10.0
- Pillow >= 9.0.0
- tkinter from the Python standard library

## Usage

### Format Conversion

1. Select a qPCR result file exported by the instrument. Roche text exports and native `.ixo` experiment files are supported.
2. Set layout parameters, including columns per sample and rows per gene.
3. Optionally enter custom sample names and gene names separated by commas.
4. Preview the parsed result or export it to Excel.
5. Use the one-click Delta Delta Ct workflow to generate side-by-side analysis output after selecting a reference gene and control group.

### Wide Table to Delta Delta Ct

Use this workflow when the input is already organized as a `Sample x Gene` wide-format Excel or CSV table.

1. Select a wide-format Excel or CSV file.
2. Select the sheet and confirm the sample/group and gene columns.
3. Choose the reference gene.
4. Select one or more target genes.
5. Choose the control group.
6. Adjust output name settings if needed.
7. Calculate Delta Delta Ct and export CSV or Excel results.

### Long-Format Delta Delta Ct

1. Select a Ct data file in Excel or CSV format.
2. Confirm column mappings.
3. Choose the primary reference gene.
4. Optionally select additional reference genes for geometric-mean normalization.
5. Choose the control group.
6. Select statistical test, figure format, and optional Excel figure embedding.
7. Calculate Delta Delta Ct and export the results.

## Output Formats

### Side-by-Side CSV

```text
Sample Name | Target Name | CT Value | Reference Mean | ... | Sample Name | Target Name | CT Value | Delta Ct | Delta Delta Ct | 2^(-Delta Delta Ct)
Control_1   | GAPDH       | 18.50    | 18.40          | ... | Control_1   | GeneA       | 24.20    | 5.70     | 0.00           | 1.00
```

### Multi-Sheet Excel

- Detailed calculation results
- Summary table
- Side-by-side analysis table
- Calculation parameters

## Formula

```text
Group reference Ct mean = mean Ct of reference genes within the same group
Sample target Ct mean = mean Ct of technical replicates for each target gene
Delta Ct = sample target Ct mean - group reference Ct mean
Control Delta Ct mean = mean Delta Ct of all control samples
Delta Delta Ct = sample Delta Ct - control Delta Ct mean
2^(-Delta Delta Ct) = relative expression or fold change
```

## Project Structure

```text
result_processor/
+-- presets/
+-- src/
|   +-- complete_gui.py
|   +-- ixo_parser.py
|   +-- plate_converter.py
+-- scripts/
|   +-- start_gui.py
+-- configs/
|   +-- example_layout.txt
|   +-- export_formats.yaml
+-- examples/
|   +-- sample_ddct_input.csv
+-- requirements.txt
+-- run_gui.bat
+-- README.md
```

## Configuration

### `configs/example_layout.txt`

```text
sample_count=8
sample_0=Control1,1,3
sample_1=Control2,4,6
...
gene_count=4
gene_0=GAPDH,A,D
gene_1=GeneA,E,H
...
```

### `configs/export_formats.yaml`

This file controls export format options, sample naming, gene naming, and UI theme settings.

## Changelog

### v2.5

- Added Pfaffl quantification support with per-gene amplification efficiency.
- Shared gene-efficiency settings across Delta Delta Ct workflows.
- Added method and efficiency records to the Excel parameter sheet.

### v2.4

- Added multi-reference gene normalization with geometric mean.
- Added statistical testing, FDR correction, and significance labels.
- Added PNG, PDF, and SVG figure export options.
- Added optional fold-change figure embedding in Excel exports.
- Added JSON analysis presets.

### v2.3

- Added reusable plate-layout presets.
- Added fold-change bar plot export.
- Added wide-table QC for missing values.
- Added automatic CSV encoding and delimiter detection.

### v2.2

- Added native `.ixo` file parsing.
- Added a wide-table-to-Delta-Delta-Ct workflow.
- Added one-click Delta Delta Ct workflow.
- Improved configuration loading and text-file encoding fallback.

### v2.1

- Added custom output filename templates.
- Added side-by-side CSV export.
- Improved Excel export with multiple sheets.

### v2.0

- Added full GUI.
- Added format conversion and Delta Delta Ct calculation.
- Added grouped analysis.

## License

MIT License
