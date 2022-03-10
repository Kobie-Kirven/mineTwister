<img src="https://github.com/Kobie-Kirven/mineTwister/blob/main/static/logo.png" width="300">
Mine a genome of interest for potential twister ribozymes

- [Installation](#installation)
    - [Install requirements](#install-requirements)
    - [Install mineTwister](#install-minetwister)
- [Tutorial](#tutorial)

## **Installation** 


### **Install requirements**

mineTwister is currently configured run on linux operating systems. This is because of it's dependency on a singularity image of R2DT. Future versions will use the Docker container of R2DT.

To run mineTwister, you need to have several dependencies installed. 

- Blast+ command line
```text
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz

tar -xzvf ncbi-blast-2.12.0+-x64-linux.tar.gz
```
Next, add the ncbi-blast-2.12.0+-x64-linux/bin folder to your path. 

- R2DT

R2DT has its own list of dependencies to run, so, we will use the singularity container instead. 

```
# Download the precompiled data
wget https://ftp.ebi.ac.uk/pub/databases/RNAcentral/r2dt/1.2/cms.tar.gz

tar -xzvf cms.tar.gz

# Build the singularity container
singularity build r2dt docker://rnacentral/r2dt
```

### **Installing mineTwister**
mineTwister can be easily installed with pip using the following command.
```
pip3 install git+https://github.com/Kobie-Kirven/mineTwister
```

And that's it! You are ready to start using mineTwister!

## **Tutorial**

In this tutorial, we will walk through a simple example using the Nipponbare IRGSP 1.0 reference genome. 

What do you need to run mineTwister?
- all you need is a genome of interest and a concensus twister sequence to search for.

We will download the reference genome using the following command:

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic.fna.gz

gunzip GCF_001433935.1_IRGSP-1.0_genomic.fna.gz
```

Great! Now we have our genome, now all we need is our consensus sequence. We can download this as well:
```
wget https://github.com/Kobie-Kirven/mineTwister/raw/main/examples/consensus.fa
```

Cool, now we have everything we need to run mineTwister. Let's make sure that mineTwister is up and running.
```
minetwister -h
usage: minetwister [-h] -i REFERENCE -c QUERY -s SINGULARITY -d DATA -o OUTPUT

Mine genomes for twister ribozymes

options:
  -h, --help      show this help message and exit
  -i REFERENCE    path to reference genome
  -c QUERY        concensus twister sequence
  -s SINGULARITY  singularity path
  -d DATA         singularity data path
  -o OUTPUT       output
```

Now let's try to run it for real. 

```
minetwister -i GCF_001433935.1_IRGSP-1.0_genomic.fna -c consensus.fa -s ./r2dt -d ./r2dt-data-cms-dev/ -o rice
```

You should see a ```rice.html``` file and when you open it, it should look like this:
<img src="https://github.com/Kobie-Kirven/mineTwister/blob/main/examples/sample_output_minetwister.png" width="300">