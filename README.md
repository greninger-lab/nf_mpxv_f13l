# NF_MPXV_F13L

MPXV consensus genome assembly and variants reporting in F13L using **ON563414.3** as reference.

## Methods

- Trim reads with fastp (default minimum length 100bp)
- Consensus genome assembly using BBMap and iVar consensus
    - minimum coverage of 10
    - minimum base quality of 20
    - minimum frequency threshold of 0.6 
- Variants reporting in F13L using iVar variants
    - minimum coverage of 10
    - minimum base quality of 20
    - minimum frequency threshold of 0.01
- Generate summary and QC stats

## Usage
Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)

Install [`Docker`](https://docs.docker.com/engine/installation/)

To run this pipeline:

	nextflow run greninger-lab/nf_mpxv_f13l -r main -latest --input example_samplesheet.csv --output example_output

with Docker:

	nextflow run greninger-lab/nf_mpxv_f13l -r main -latest --input example_samplesheet.csv --output example_output -profile docker

on AWS:
    
	nextflow run greninger-lab/nf_mpxv_f13l -r main -latest --input example_samplesheet.csv --output example_output -profile docker -c your_nextflow_aws.config
	

## Options
|Option|Explanation|
|------|-----------|
| `--input` | samplesheet in csv format with fastq information |
| `--output` | output directory (default: nf_mpxv_f13l_output) |
| `--run_name` | name for the summary tsv file (default: 'run') |
| `--skip_fastqc` | skip quality control using FastQC (default: false) |
| `--skip_fastp` | skip adapters and reads trimming using fastp (default: false) |
| `--trim_len` | minimum read length to keep (default:100) |
| `--save_trimmed_reads` | save trimmed fastq (default: false) |
| `--sample` | downsample fastq to a certain fraction or number of reads |
| `--min_mapped_reads` | minimum number of mapped reads for variants and consensus calling (default: 1000) |
| `--ivar_variants_t` | minimum frequency threshold to call consensus (default: 0.01) |
| `--ivar_variants_q` | minimum quality score threshold to call consensus (default: 20) |
| `--ivar_variants_m` | minimum depth to call consensus (default: 10) |
| `--save_mpileup` | save samtools mpileup used for iVar variants |
| `--ivar_consensus_t` | minimum frequency threshold to call consensus (default: 0.6) |
| `--ivar_consensus_q` | minimum quality score threshold to call consensus (default: 20) |
| `--ivar_consensus_m` | minimum depth to call consensus (default: 10) |

## Usage notes
- Samplesheet example: `assets/samplesheet.csv`
- You can create a samplesheet using the bundled python script: `python bin/fastq_dir_samplesheet.py fastq_dir samplesheet_name.csv`
- Memory and CPU usage for pipeline processes can be adjusted in `conf/base.config`
- Process arguments can be adjusted in `conf/modules.config`
- If you are using Docker on Linux, check out these [post-installation steps](https://docs.docker.com/engine/install/linux-postinstall/) (especially cgroup swap limit capabilities support) for configuring Linux to work better with Docker.- By default, Docker has full access to full RAM and CPU resources of the host, but if you are using MacOS, go to Settings -> Resources in Docker Desktop to make sure enough resources are allocated to docker containers. 

## Contact
For bug reports please email aseree@uw.edu or raise an issue on Github.
