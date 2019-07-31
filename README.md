<a href="https://www.sentieon.com/">		<img src="https://www.sentieon.com/wp-content/uploads/2017/05/cropped-companylogo.png"  alt="Sentieon" width="25%" >	</a>

# Run Sentieon pipelines on Google Cloud Platform

Here we present a convenient tool for your to run our typical pieplines on Google Cloud Plarform. If you want to customize your pipeline, you may want to visit https://www.sentieon.com/support/.  

## Pipelines

- Mapping(FASTQ -> BAM): `sentieon bwa-mem`
- Variant Calling(BAM -> VCF) 
  - Germline: `DNAseq` `DNAscope`
    - WGS
    - WES
  - Somatic: `TNhaplotyper` `TNhaplotyper2` `TNscope`
    - Tumor-Normal
    - Tumor only
    
## Prerequisite & Setup

- Python 2.7+
- GCP account and Cloud SDK

  You may want to follow the instructions in "Before you begin" section in this [Google Cloud Genomics Tutorial](https://cloud.google.com/genomics/docs/tutorials/sentieon) to set up the working environment.
- Set up your virtual environment and install dependencies
   ```bash
  #pip install virtualenv
  virtualenv env
  source env/bin/activate
  pip install --upgrade \
      pyyaml \
      google-api-python-client \
      google-auth \
      google-cloud-storage \
      google-auth-httplib2
   ```
- Download the example files and set your current directory

   ```bash
   git clone https://github.com/sentieon/sentieon-google-genomics.git
   cd sentieon-google-genomics
   ```

## Run your first Sentieon job on GCP

Right now, we are granting free-trial license to your account automatically. You will get 14 days free trial starting from the time you run your first Sentieon job.

  1. Specify your `PROJECT` and `BUCKET` in `examples/example.json`
  ```json
  {
  "FQ1": "gs://sentieon-test/pipeline_test/inputs/test1_1.fastq.gz",
  "FQ2": "gs://sentieon-test/pipeline_test/inputs/test1_2.fastq.gz",
  "REF": "gs://sentieon-test/pipeline_test/reference/hs37d5.fa",
  "OUTPUT_BUCKET": "YOUR_BUCKET_HERE",
  "ZONES": "us-central1-a,us-central1-b,us-central1-c,us-central1-f",
  "PROJECT_ID": "YOUR_PROJECT_HERE"
  }
  ```

  2. Run
  ```bash
  python runner/sentieon_runner.py examples/example.json
  ```

## Pipeline Examples

   ### JSON examples

   In `examples` directory, you could find 
   
   file | pipeline 
   --- | ---
   100x_wes.json | DNAseq pipeline from FASTQ to VCF for Whole Exome Sequencing Data
   30x_wgs.json | DNAseq pipeline from FASTQ to VCF for Whole Genome Sequencing Data
   tn_example.json | TNseq pipeline from FASTQ to VCF for Tumor Normal Pairs
   
   ### Parameters
   
   In `runner/sentieon_germline.yaml` and `runner/sentieon_somatic.yaml`, you could find adjustable parameters for germline and somatic pipelines.
   
   

  
