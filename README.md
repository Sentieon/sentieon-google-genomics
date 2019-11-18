<a href="https://www.sentieon.com/">  <img src="https://www.sentieon.com/wp-content/uploads/2017/05/cropped-companylogo.png" alt="Sentieon" width="25%">	</a>

# Run Sentieon pipelines on the Google Cloud Platform
*Easy to use pipelines for the Sentieon tools on the Google Cloud*

For a tutorial, see Google's tutorial on [running a Sentieon DNAseq pipeline](https://cloud.google.com/genomics/docs/tutorials/sentieon). For more customized pipelines and additional details on the Sentieon software, please visit https://www.sentieon.com.

## Table of Contents
- [Highlights](#highlights)
- [Prerequisites](#prerequisites)
- [Running a pipeline](#running)
  - [Configure your environment](#configure)
  - [Download the example files](#example)
  - [Understanding the input format](#input)
  - [Run the example pipeline](#run)
  - [Understanding the output](#understand)
  - [Other example pipelines](#other_examples)
- [Recommended configurations](#recommended)
  - [Germline WGS and WES](#germline)
  - [Somatic WGS and WES](#somatic)
- [Additional options - germline](#configurations_germline)
  - [Input file options](#germline_input)
  - [Machine options](#germline_machine)
  - [Pipeline configuration](#germline_config)
  - [Pipeline options](#germline_options)
- [Additional options - somatic](#configurations_somatic)
  - [Input file options](#somatic_input)
  - [Machine options](#somatic_machine)
  - [Pipeline configuration](#somatic_config)
  - [Pipeline options](#somatic_options)
- [Where to get help](#help)

<a name="highlights"/>

## Highlights

- Easily run the Sentieon Pipelines on the Google Cloud.
- Pipelines are optimized by Sentieon to be well-tuned for efficiently processing WES and WGS samples.
- All Sentieon pipelines and variant callers are available including DNAseq, DNAscope, TNseq, and TNscope
- Matching results to the GATK [Germline](https://github.com/Sentieon/sentieon-dnaseq) and Somatic Best Practices Pipelines
- Automatic 14-day free-trial of the Sentieon software on the Google Cloud

<a name="prerequisites"/>

## Prerequisites

1. [Install Python 2.7+](https://www.python.org/downloads/).
2. Select or create a GCP project.
3. Make sure that billing is enabled for your Google Cloud Platform project.
4. Enable the Cloud Life Sciences, Compute Engine, and Cloud Storage APIs.
5. Install and initialize the Cloud SDK.
6. Update and install gcloud components:
```bash
gcloud components update &&
gcloud components install alpha
```
7. [Install git](https://git-scm.com/downloads) to download the required files.
8. By default, Compute Engine has resource quotas in place to prevent inadvertent usage. By increasing quotas, you can launch more virtual machines concurrently, increasing throughput and reducing turnaround time.
For best results in this tutorial, you should request additional quota above your project's default. Recommendations for quota increases are provided in the following list, as well as the minimum quotas needed to run the tutorial. Make your quota requests in the us-central1 region:
  CPUs: 64
  Persistent Disk Standard (GB): 375
You can leave other quota request fields empty to keep your current quotas.

<a name="running"/>

## Running a pipeline


<a name="configure"/>

### Configuring your environment

Setup a Python virtualenv to manage the environment. First, install virtualenv if necessary
```bash
pip install --upgrade virtualenv
```
Install the required Python dependencies
```bash
virtualenv env
source env/bin/activate
pip install --upgrade \
    pyyaml \
    google-api-python-client \
    google-auth \
    google-cloud-storage \
    google-auth-httplib2
```

<a name="example"/>

### Download the pipeline script

Download the pipeline script and move into the new directory.
```bash
git clone https://github.com/sentieon/sentieon-google-genomics.git
cd sentieon-google-genomics
```

<a name="input"/>

### Understanding the input format

The runner script accepts a JSON file as input. In the repository you downloaded, there is an `examples/example.json` file with the following content:
```json
{
  "FQ1": "gs://sentieon-test/pipeline_test/inputs/test1_1.fastq.gz",
  "FQ2": "gs://sentieon-test/pipeline_test/inputs/test1_2.fastq.gz",
  "REF": "gs://sentieon-test/pipeline_test/reference/hs37d5.fa",
  "OUTPUT_BUCKET": "gs://BUCKET",
  "ZONES": "us-central1-a,us-central1-b,us-central1-c,us-central1-f",
  "PROJECT_ID": "PROJECT_ID",
  "EMAIL": "EMAIL"
}
```

The following table describes the JSON keys in the file:

| JSON key      | Description                                                                   |
| ------------- | ----------------------------------------------------------------------------- |
| FQ1           | The first pair of reads in the input fastq file.                              |
| FQ2           | The second pair of reads in the input fastq file.                             |
| BAM           | The input BAM file, if applicable.                                            |
| REF           | The reference genome. If set, the reference index files are assumed to exist. |
| OUTPUT_BUCKET | The bucket and directory used to store the data output from the pipeline.     |
| ZONES         | A comma-separated list of GCP zones to use for the worker node.               |
| PROJECT_ID    | Your GCP project ID.                                                          |
| EMAIL         | Your email                                                                    |

The `FQ1`, `FQ2`, `REF`, and `ZONES` fields will work with the defaults. However, the `OUTPUT_BUCKET`, `PROJECT_ID`, and `EMAIL` fields will need to be updated to point to your specific output bucket/path, Project ID, and email address.

<a name="run"/>

### Run the example pipelines

Edit the `OUTPUT_BUCKET`, `PROJECT_ID`, and `EMAIL` fields in the `examples/example.json` to your output bucket/path, the GCP Project ID that you setup earlier, and email you want associated with your Sentieon license. By supplying the `EMAIL` field, your PROJECT_ID will automatically receive a 14 day free trial for the Sentieon software on the Google Cloud.

You after modifying the `examples/example.json` file, you can use the following command to run the DNAseq pipeline on a small test dataset.
```bash
python runner/sentieon_runner.py examples/example.json
```

<a name="understand"/>

### Understanding the output

If execution is successful, the runner script will print some logging information followed by `Operation succeeded` to the terminal. Output files from the pipeline can then be found in the `OUTPUT_BUCKET` location in Google Cloud Storage including alignment (BAM) files, variant calls, sample metrics and logging information.

In the event of run failure, some diagnostic information will be printed to the screen followed by an error message. For assistance, please send the diagnostic information along with any log files in `OUTPUT_BUCKET`/worker_logs/ to support@sentieon.com.

<a name="other_examples"/>

### Other example pipelines
In the `examples` directory, you can find the following example configurations:

| Configuration   | Pipeline                                                           |
| --------------- | ------------------------------------------------------------------ |
| 100x_wes.json   | DNAseq pipeline from FASTQ to VCF for Whole Exome Sequencing Data  |
| 30x_wgs.json    | DNAseq pipeline from FASTQ to VCF for Whole Genome Sequencing Data |
| tn_example.json | TNseq pipeline from FASTQ to VCF for Tumor Normal Pairs            |

<a name="recommended"/>

## Recommended configurations

Below are some recommended configurations for some common use-cases. The cost and runtime estimates below assume that jobs are run on preemptible instances that were not preempted during job execution.

<a name="germline"/>

### Germline Whole Genome and Whole Exome Sequencing

The following configuration will run a 30x human genome at a cost of approximately $1.35 and will take about 2 hours. This configuration can also be used to run a 100x whole exome at a cost of approximately $0.30 and will take about 35 minutes.
```json
{
  "FQ1": "gs://my-bucket/sample1_1.fastq.gz",
  "FQ2": "gs://my-bucket/sample1._2.fastq.gz",
  "REF": "gs://sentieon-test/pipeline_test/reference/hs37d5.fa",
  "OUTPUT_BUCKET": "gs://BUCKET",
  "ZONES": "us-central1-a,us-central1-b,us-central1-c,us-central1-f",
  "PROJECT_ID": "PROJECT_ID",
  "EMAIL": "EMAIL",
  "BQSR_SITES": "gs://sentieon-test/pipeline_test/reference/Mills_and_1000G_gold_standard.indels.b37.vcf.gz,gs://sentieon-test/pipeline_test/reference/1000G_phase1.indels.b37.vcf.gz,gs://sentieon-test/pipeline_test/reference/dbsnp_138.b37.vcf.gz",
  "DBSNP": "gs://sentieon-test/pipeline_test/reference/dbsnp_138.b37.vcf.gz",
  "PREEMPTIBLE_TRIES": "2",
  "NONPREEMPTIBLE_TRY": true,
  "STREAM_INPUT": "True",
  "DISK_SIZE": 300,
  "PIPELINE": "GERMLINE",
  "CALLING_ALGO": "Haplotyper"
}
```
The `CALLING_ALGO` key can be changed to `DNAscope` to use the Sentieon DNAscope variant caller for improved variant calling accuracy. For large input files, `DISK_SIZE` should be increased.


<a name="somatic"/>

### Somatic Whole Genome and Whole Exome Sequencing

The following configuration will run a paired 60-30x human genome at a cost of approximately $3.70 and will take about 7 hours. This configuration can also be used to run a paired 150-150x human exome at a cost of approximately $0.60 and will take about 1.5 hours.
```json
{
  "TUMOR_FQ1": "gs://my-bucket/tumor1_1.fastq.gz",
  "TUMOR_FQ2": "gs://my-bucket/tumor1_2.fastq.gz",
  "FQ1": "gs://my-bucket/normal1_1.fastq.gz",
  "FQ2": "gs://my-bucket/normal1._2.fastq.gz",
  "REF": "gs://sentieon-test/pipeline_test/reference/hs37d5.fa",
  "OUTPUT_BUCKET": "gs://BUCKET",
  "ZONES": "us-central1-a,us-central1-b,us-central1-c,us-central1-f",
  "PROJECT_ID": "PROJECT_ID",
  "EMAIL": "EMAIL",
  "BQSR_SITES": "gs://sentieon-test/pipeline_test/reference/Mills_and_1000G_gold_standard.indels.b37.vcf.gz,gs://sentieon-test/pipeline_test/reference/1000G_phase1.indels.b37.vcf.gz,gs://sentieon-test/pipeline_test/reference/dbsnp_138.b37.vcf.gz",
  "DBSNP": "gs://sentieon-test/pipeline_test/reference/dbsnp_138.b37.vcf.gz",
  "PREEMPTIBLE_TRIES": "2",
  "NONPREEMPTIBLE_TRY": true,
  "STREAM_INPUT": "True",
  "DISK_SIZE": 300,
  "PIPELINE": "SOMATIC",
  "CALLING_ALGO": "TNhaplotyper"
}
```
The `CALLING_ALGO` key key can be change to `TNsnv`, `TNhaplotyper`, `TNhaplotyper2`, or `TNscope` to use Sentieon's TNsnv, TNhaplotyper, TNhaplotyper2 or TNscope variant callers, respectively. For large input files, `DISK_SIZE` should be increased.

<a name="configurations_germline"/>

## Additional options - Germline

<a name="germline_input"/>

### Input file options

| JSON Key       | Description                                                                          |
|--------------- | ------------------------------------------------------------------------------------ |
| FQ1            | A comma-separated list of input R1 FASTQ files                                       |
| FQ2            | A comma-separated list of input R2 FASTQ files                                       |
| READGROUP      | A comma-separted list of readgroups headers to add to the read data during alignment |
| BAM            | A comma-separated list of input BAM files                                            |
| REF            | The path to the reference genome                                                     |
| BQSR_SITES     | A comma-separated list of known sites for BQSR                                       |
| DBSNP          | A dbSNP file to use during variant calling                                           |
| INTERVAL       | A string of interval(s) to use during variant calling                                |
| INTERVAL_FILE  | A file of intervals(s) to use during variant calling                                 |
| DNASCOPE_MODEL | A trained model to use during DNAscope variant calling                               |

<a name="germline_machine"/>

### Machine options

| JSON Key     | Description                                                                 |
| ------------ | --------------------------------------------------------------------------- |
| ZONES        | GCE Zones to potentially launch the job in                                  |
| DISK_SIZE    | The size of the hard disk to use (should be 3x the size of the input files) |
| MACHINE_TYPE | The type of GCE machine to use to run the pipeline                          |

<a name="germline_config"/>

### Pipeline configurations

| JSON Key            | Description                                                             |
| ------------------- | ----------------------------------------------------------------------- |
| SENTIEON_VERSION    | The version of the Sentieon software package to use                     |
| DEDUP               | Type of duplicate removal to run (nodup, markdup or rmdup)              |
| NO_METRICS          | Skip running metrics collection                                         |
| NO_BAM_OUTPUT       | Skip outputting a preprocessed BAM file                                 |
| NO_HAPLOTYPER       | Skip variant calling                                                    |
| GVCF_OUTPUT         | Output variant calls in gVCF format rather than VCF format              |
| STREAM_INPUT        | Stream the input FASTQ files directly from Google Cloud Storage         |
| RECALIBRATED_OUTPUT | Apply BQSR to the output preprocessed alignments (not recommended)      |
| CALLING_ARGS        | A string of additional arguments to pass to the variant caller          |
| PIPELINE            | Set to `GERMLINE` to run the germline variant calling pipeline          |
| CALLING_ALGO        | The Sentieon variant calling algo to run. Either Haplotyper or DNAscope |

<a name="germline_options"/>

### Pipeline options

| JSON Key            | Description                                                                                         |
| ------------------- | --------------------------------------------------------------------------------------------------- |
| OUTPUT_BUCKET       | The Google Cloud Storage Bucket and path prefix to use for the output files                         |
| EMAIL               | An email address to use to obtain an evaluation license for your GCP Project                        |
| SENTIEON_KEY        | Your Sentieon license key (only applicable for paying customers)                                    |
| PROJECT_ID          | Your GCP Project ID to use when running jobs                                                        |
| PREEMPTIBLE_TRIES   | Number of attempts to run the pipeline using preemptible instances                                  |
| NONPREEMPTIBLE_TRY  | After `PREEMPTIBLE_TRIES` are exhausted, whether to try one additional run with standard instances  |

<a name="configurations_somatic"/>

## Additional options - Somatic

<a name="somatic_input"/>

### Input file options

| JSON Key        | Description                                                                                 |
|---------------- | ------------------------------------------------------------------------------------------- |
| TUMOR_FQ1       | A comma-separated list of input R1 tumor FASTQ files                                        |
| TUMOR FQ2       | A comma-separated list of input R2 tumor FASTQ files                                        |
| FQ1             | A comma-separated list of input R1 normal FASTQ files                                       |
| FQ2             | A comma-separated list of input R2 normal FASTQ files                                       |
| TUMOR_READGROUP | A comma-separted list of readgroups headers to add to the tumor read data during alignment  |
| READGROUP       | A comma-separted list of readgroups headers to add to the normal read data during alignment |
| TUMOR_BAM       | A comma-separated list of input tumor BAM files                                             |
| BAM             | A comma-separated list of input normal BAM files                                            |
| REF             | The path to the reference genome                                                            |
| BQSR_SITES      | A comma-separated list of known sites for BQSR                                              |
| DBSNP           | A dbSNP file to use during variant calling                                                  |
| INTERVAL        | A string of interval(s) to use during variant calling                                       |
| INTERVAL_FILE   | A file of intervals(s) to use during variant calling                                        |

<a name="somatic_machine"/>

### Machine options

| JSON Key     | Description                                                                 |
| ------------ | --------------------------------------------------------------------------- |
| ZONES        | GCE Zones to potentially launch the job in                                  |
| DISK_SIZE    | The size of the hard disk to use (should be 3x the size of the input files) |
| MACHINE_TYPE | The type of GCE machine to use to run the pipeline                          |

<a name="somatic_config"/>

### Pipeline configurations

| JSON Key            | Description                                                                                             |
| ------------------- | ------------------------------------------------------------------------------------------------------- |
| SENTIEON_VERSION    | The version of the Sentieon software package to use                                                     |
| DEDUP               | Type of duplicate removal to run (nodup, markdup or rmdup)                                              |
| NO_METRICS          | Skip running metrics collection                                                                         |
| NO_BAM_OUTPUT       | Skip outputting a preprocessed BAM file                                                                 |
| NO_VCF              | Skip variant calling                                                                                    |
| STREAM_INPUT        | Stream the input FASTQ files directly from Google Cloud Storage                                         |
| RECALIBRATED_OUTPUT | Apply BQSR to the output preprocessed alignments (not recommended)                                      |
| CALLING_ARGS        | A string of additional arguments to pass to the variant caller                                          |
| PIPELINE            | Set to `SOMATIC` to run the somatic variant calling pipeline                                            |
| RUN_TNSNV           | If using the TNseq pipeline, use TNsnv for variant calling                                              |
| CALLING_ALGO        | The Sentieon somatic variant calling algo to run. Either TNsnv, TNhaplotyper, TNhaplotyper2, or TNscope |

<a name="somatic_options"/>

### Pipeline options

| JSON Key            | Description                                                                                         |
| ------------------- | --------------------------------------------------------------------------------------------------- |
| OUTPUT_BUCKET       | The Google Cloud Storage Bucket and path prefix to use for the output files                         |
| EMAIL               | An email address to use to obtain an evaluation license for your GCP Project                        |
| SENTIEON_KEY        | Your Sentieon license key (only applicable for paying customers)                                    |
| PROJECT_ID          | Your GCP Project ID to use when running jobs                                                        |
| PREEMPTIBLE_TRIES   | Number of attempts to run the pipeline using preemptible instances                                  |
| NONPREEMPTIBLE_TRY  | After `PREEMPTIBLE_TRIES` are exhausted, whether to try one additional run with standard instances  |

<a name="help"/>

## Where to get help

Please email support@sentieon.com with any questions.
