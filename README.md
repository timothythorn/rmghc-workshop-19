# Rocky Mountain Genomics Hackcon Workshop 2019

An introduction to tools for reproducible analysis pipelines.

## Prerequisites
- a text editor ([Atom](https://atom.io/) recommended)
- sftp client ([FileZilla](https://filezilla-project.org/download.php) recommended)
- a terminal/ssh client ([Putty](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html) for windows; terminal or [iTerm](https://www.iterm2.com/) for mac)
- a GitHub account

## Connect to your Google Cloud Compute Instance

**Challenge round:** What operating system is your compute instance running?

## Set up your git repository

Fork this repository so you have your own copy

Clone your new repository
```bash
git clone https://github.com/msmallegan/rmghc-workshop-19.git
```


## Edit files two ways

Use your file manager to edit a file on your compute instance on your local machine.

Use nano on the command line to edit the same file. 

Commit your changes and push to GitHub

## Interact with singularity

Let's inspect the singularity container. What software is loaded on it?

*Run something on singularity with singularity exec*

## Run the RNA-seq analysis pipeline with nextflow

```bash
nextflow run main.nf \
  -resume \
  -with-report ./reports/rnaseq_test.html \
  -with-dag ./reports/rnaseq_test_dag.pdf
```

## View the pipeline reports in your browser

## Add a process to run the differential expression RScript

## View the results of the RMarkdown script
