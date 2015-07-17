#DNMFilter
#Version 0.1.1

1. Installation 

1) Install Java JDK or JRE; 

2) Install R (Rscript exectuable must be on the path); 

3) Install Runiversal package(http://cran.rproject.org/web/packages/Runiversal/index.html) in R environment; 

4) Install gbm package (http://cran.rproject.org/web/packages/gbm/index.html) in R environment.

2. Running 

usage: java -jar DNMFilter.jar <COMMAND> [OPTIONS] 

1) extract 

This command extracts sequence features of known true or false positive DNMs to build the training set: 

usage: java -jar DNMFilter.jar extract [OPTIONS] 

--reference	<FILE> reference genome file (required) 
--pedigree		<FILE> pedigree file (required) 
--positive		<FILE> known true positive DNM file (required) 
--negative		<FILE> known false positive DNM file (required) 
--bam			<FILE> bam list file (required) 
--output		<FILE> output file (required) 

2) gbm 

This command uses gradient boosting approach to filter de novo mutations: 

usage: java -jar DNMFilter.jar gbm [OPTIONS] 

--reference	<FILE> reference genome file (required) 
--pedigree		<FILE> pedigree file (required) 
--bam			<FILE> bam list file (required) 
--training		<FILE> training set (required) 
--candidate	<FILE> candidate DNM file (required) 
--configuration	<FILE> feature configuration file (required) 
--cutoff		<DOUBLE> cutoff to determine a putative DNM (optional, default 0.4) 
--output		<FILE> output file (required) 

3. File Instruction

1) bam list file (two columns, tab-separated) 

Column 1: sampleID 
Column 2: path of .bam file

Example: 

Sample1	/path/Sample1.bam 
Sample2	/path/Sample2.bam 
Sample3	/path/Sample3.bam 
Sample4	/path/Sample4.bam
...		...	

2) DNM file (.csv file, the first three columns must be specified, comma-separated) 

Column 1: familyID 
Column 2: chromsome 
Column 3: position

Example: 

Trio1,1,10000 
Trio1,2,20000 
Trio1,2,30000 
Trio2,1,40000

Note: The DNM file should be first sorted by familyID and then by chromosome position.

3) feature configuration file

The feature that values 1 is selected to train the model and filter DNMs, while the feature that values 0 is not selected.

4) pedigree file

See (http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml)

Note the Family ID must be the same as the first column of DNM file and the Individual ID must be the same as the first column of bam list file

4. Others

The package includes a training set built with 264 Epi4k exome trios, which users can employ it to filter new potential DNMs. In addition, the package also includes a feature configuration file.

5. Contact
yzhuangliu@gmail.com;  yongzhuang.liu@hit.edu.cn
