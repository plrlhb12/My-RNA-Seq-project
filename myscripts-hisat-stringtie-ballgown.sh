# generate a enviroment which contains a new folder $Home/bin
# later I will copy my exectuble files inside this enviroment therefore the computer will search it auto
$ mkdir $HOME/bin
$ expot PATH=$HOME/bin:$PATH

# cd to the workind directory, where all files of data and analysis are in it
$ cd ~/Downloads/Bioinformatics-RNA-seq/ftp.ccb.jhu.edu/pub/RNASeq_protocal
# get data files which also contain genenames, geneids, indexes and samples
$ wget -c -r ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol
# unzip the datafiles
$ tar xvzf chrX_data.tar.gz
$ brew install tree
$ tree chrX_data


# pre-download the package of samtools-1.10 into the working directory and 
$ tar jxvf samtools-1.10.tar.bz2
# find the exeutable file samtools-1
$ cd samtools-1.10
$ ./configure --prefix=$HOME/bin
$ make
$ make install
# copy it to the enviroment $HOME/bin
$ cp samtools-1.10/samtools-1 $HOME/bin

# download the package of hisat2
$ wget wget -c ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-OSX_x86_64.zip
$ unzip hisat2-2.1.0-OSX_x86_64.zip
$ cp hisat2* *.py $HOME/bin

$ wget -c http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3b.OSX_x86_64.tar
$ tar xvzf stringtie-2.1.0.OSX_x86_64.tar.gz 
$ cp stringtie-2.1.0.OSX_x86_64/stringtie $HOME/bin

# generate shell file for batch executation
$ vi peng_map.sh   # use hisat2
$ sh peng_map.sh
$ vi peng_convert_bam.sh  # use samtool
$ sh peng_convert_bam.sh
$ vi peng_assemble.sh   # use stringtie
$ sh peng_assemble.sh

# do merage
$ stringtie --merge -p 8 -G chrX_data/genes/chrX.gtf -o stringtie_merged.gtf chrX_data/mergelist.txt
$ cat stringtie_merged.gtf | grep -v "^#" | awk '$3=="transcript" {print}' | wc -l

# compared to annotated, using gffcompare
$ git clone https://github.com/gpertea/gclib
$ git clone https://github.com/gpertea/gffcompare
$ cd gffcompare
$ make release
$ cp gffcompare $HOME/bin
$ gffcompare -r chrX_data/genes/chrX.gtf -G -o merged stringtie_merged.gtf
$ cat merged.stats

$ vi peng_ballgown.sh
$ sh peng_ballgown.sh

# check other files in chrx_data
$ cd chrx_data
$ less geuvadis_phenodata.csv
$ less mergelist.txt










