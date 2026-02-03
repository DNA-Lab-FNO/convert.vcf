# convert.vcf

Conversion of VCF to FinalistDX Excel format. This package depends on a
[custom version of vcf.R package](https://github.com/gorgitko/vcfR) that can parse annotated VCFs produced by
VEP and snpEff tools.

## Installation

### Package

```r
install.packages("remotes")
# Replace "<PAT>" with your GitHub PAT
Sys.setenv(GITHUB_PAT = "<PAT>")
remotes::install_github("DNA-Lab-FNO/convert.vcf")
```

### Command line interface (CLI)

This packages provides a simple CLI script that calls `convert_vcf_files_to_finalist()` followed by `write_output_file()`.
To install the script on your system, use

```r
convert.vcf::install_cli("/path/to/directory")
```

For example, you can use:

- `/usr/local/bin` for system-wide installation
- `~/.local/bin` for user-specific installation
  - Note that this directory must be present in your `PATH` environment variable.
    If not, run `export PATH=${PATH}:~/.local/bin`, or put this line in your `~/.bashrc` to make it permanent.
- `/opt/miniconda3/envs/<env>/bin` for installation into a specific conda environment

## Usage

### In R

After installation see help pages of the functions below (e.g. call `?read_vcf` in R)

- `read_vcf()`
- `convert_vcf_files_to_finalist()`

### From CLI

For available arguments, run `convert_vcf.R -h`

Example of batch conversion of VCF files into FinalistDX Excel format
(note that `haplotypecaller` and `vep` are default expected variant calling and variant annotation tools, respectively):

```bash
convert_vcf.R \
  --output_file 2024_07_03_Panel_run_4_BRCA202.haplotypecaller.filtered_VEP.xlsx \
  /data/scratch/run_sarek/2024_12_04/2024_07_03_Panel_run_4/BRCA202/results/variant_calling/results/annotation/haplotypecaller/*/*.haplotypecaller.filtered_VEP.ann.vcf.gz
```
