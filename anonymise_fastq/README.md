----
# Anonymise FastQ files
----
We wish to strip away any personal information contained in a FastQ ID entry (the lines with "@") while retaining information for downstream bioinformatic analysis. 

This module will provide the facility to stream FastQ entries from (perhaps sense-antisense pairs of) FastQ files and output anonymised files.  

There are few things to be decided : 
  - What are the information that need to persist? 
    - Perhaps just "lane" information. 
  - How do we parse the original ID information? 
    - The ID line format are sequencing provider dependent. 
    - We will provide a highly costumisable anonymising function. 
  - What anonymisation procedure do we use? 


### The anonymising function
----
Using a specified id format, we will first build a dictionary `{field_name:field_value}` pairs that has all the information contained in the input id.  

Then given a list of desirable fields to be retained, we will construct a new anonymised id where the desired information can still be recoverable. 



### Tracing Memory Footprint
---
We will use Python `memory_profiling` package to trace memory usage to make sure that we are not reading entire FastQ file into memory. 





----
# Install
----
1. Run `$ git clone https://github.com/edmundlth/MGHA_bioinformatics`. 
2. Run `$ pip install MGHA_bioinformatics/anonymise_fastq`. 
3. See `$ anonymise_fastq --help` for input options. 

----
# Example run
----
### On commandline
----
1. Using default settings (need to make sure ID format match default. See `anonymise_fastq --help`): 
`$ anonymise_fastq --infilepath /path/to/input.fastq.gz`
2. Displaying all options: 
`$ anonymise_fastq --infilepath testinput.fastq.gz --outfilepath testoutput.fastq.gz --idformat <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>/<sense> --retained_fields lane title`

### Single Input
----
We will be processing a single FastQ file with id format given by  

```<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>/<sense>```  

e.g. 

```DAJHDKJHF78Q1:270:H291LADAX:1:1206:8013:75391/1```. 

For demo purpose, we will retain the following fields: 

`"lane", "title","sense"`

separated by `:` and generate a (probably random) new name by prefixing the resulting string with current unix time. With the sample id above, the output will be

`'1554777537:1:1206:1'`. 