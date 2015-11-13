#!/usr/bin/perl
use strict;
use warnings;
use Compress::Zlib;
#####################
# Kaanan Shah
# 11/13/15
# script to create glmnet predictors of gene expression from GTEx data
# This script sets up multiple jobs on tarbell - 1 per chr per tissue
# The second step of the script will concatenate all the per chr files into one file

my $dotissuewide = 0;
my $cattissuewide = 0;

my @tissues = ("Adipose-Subcutaneous","AdrenalGland","Artery-Aorta","Artery-Coronary","Artery-Tibial",
"Brain-Anteriorcingulatecortex(BA24)","Brain-Caudate(basalganglia)","Brain-CerebellarHemisphere","Brain-Cerebellum",
"Brain-Cortex","Brain-FrontalCortex(BA9)","Brain-Hippocampus","Brain-Hypothalamus","Brain-Nucleusaccumbens(basalganglia)",
"Brain-Putamen(basalganglia)","Breast-MammaryTissue","Cells-EBV-transformedlymphocytes","Cells-Transformedfibroblasts",
"Colon-Sigmoid","Colon-Transverse","Esophagus-GastroesophagealJunction","Esophagus-Mucosa","Esophagus-Muscularis","Liver",
"Lung","Muscle-Skeletal","Nerve-Tibial","Ovary","Pancreas","Pituitary","Skin-NotSunExposed(Suprapubic)","Skin-SunExposed(Lowerleg)","SmallIntestine-TerminalIleum","Spleen","Stomach","Testis","Thyroid","WholeBlood","Heart-AtrialAppendage","Heart-LeftVentricle");


my @chromosomes = 1 .. 22;
my @alphas = (0.5); #0,0.05,0.95,1,0.5);

if ($dotissuewide == 1) {
    foreach my $a (@alphas) {
        foreach my $t (@tissues) {
            foreach my $chrom (@chromosomes) {
                ## deal with the annoying parentheses in the file names... who does that?! come on man! unix hates special characters!
                my $filet = $t;
                if ($filet eq "Brain-Anteriorcingulatecortex(BA24)") {$filet = "Brain-Anteriorcingulatecortex-BA24";}
                if ($filet eq "Brain-Caudate(basalganglia)") {$filet = "Brain-Caudate-basalganglia";}
                if ($filet eq "Brain-FrontalCortex(BA9)") {$filet = "Brain-FrontalCortex-BA9";}
                if ($filet eq "Brain-Nucleusaccumbens(basalganglia)") {$filet = "Brain-Nucleusaccumbens-basalganglia";}
                if ($filet eq "Brain-Putamen(basalganglia)") {$filet = "Brain-Putamen-basalganglia";}
                if ($filet eq "Skin-NotSunExposed(Suprapubic)") {$filet = "Skin-NotSunExposed-Suprapubic";}
                if ($filet eq "Skin-SunExposed(Lowerleg)") {$filet = "Skin-SunExposed-Lowerleg";}
                
                open (R, ">/group/im-lab/nas40t2/kaanan/PrediXcan/betas/scripts/TW_${filet}_${chrom}_${a}.R") or die "cant make TW_${filet}_${chrom}_${a}.R";
                print R "source(\"/group/im-lab/nas40t2/kaanan/PrediXcan/betas/scripts/26_GTEx_TW_CV_elasticNet_function.r\")\n";
                print R "doTW(\"${t}\",${chrom},${a})\n";
                close(R);
                open (RUN, ">/group/im-lab/nas40t2/kaanan/PrediXcan/betas/scripts/TW_${filet}_${chrom}_${a}.txt") or die "cant make /group/im-lab/nas40t2/kaanan/PrediXcan/betas/scripts/TW_${filet}_${chrom}_${a}.txt\n";
                print RUN "#!/bin/bash\n";
                print RUN "#PBS -N ${filet}_${chrom}_${a}_TW\n";
                print RUN "#PBS -S /bin/bash\n";
                print RUN "#PBS -l mem=5gb\n";
                print RUN "#PBS -o /group/im-lab/nas40t2/kaanan/PrediXcan/betas/scripts/logs/TW_${filet}_${chrom}_${a}.out\n";
                print RUN "#PBS -e /group/im-lab/nas40t2/kaanan/PrediXcan/betas/scripts/logs/TW_${filet}_${chrom}_${a}.err\n";
                print RUN "#PBS -l walltime=100:00:00\n";
                print RUN "module load R/3.1.0\n";
                print RUN "R --no-save < /group/im-lab/nas40t2/kaanan/PrediXcan/betas/scripts/TW_${filet}_${chrom}_${a}.R\n";
                close(RUN);
                system("qsub /group/im-lab/nas40t2/kaanan/PrediXcan/betas/scripts/TW_${filet}_${chrom}_${a}.txt");
            }
        }
    }
}

if ($cattissuewide == 1) {
    foreach my $a (@alphas) {
        system("mkdir /group/im-lab/nas40t2/kaanan/PrediXcan/betas/allBetas/elasticNet${a}/");
        system("mkdir /group/im-lab/nas40t2/kaanan/PrediXcan/betas/allBetas/PrediXcanFormat/elasticNet${a}/");
        system("mkdir /group/im-lab/nas40t2/kaanan/PrediXcan/betas/allResults/elasticNet${a}/");
        foreach my $t (@tissues) {
            my $filet = $t;
            if ($filet eq "Brain-Anteriorcingulatecortex(BA24)") {$filet = "Brain-Anteriorcingulatecortex-BA24";}
            if ($filet eq "Brain-Caudate(basalganglia)") {$filet = "Brain-Caudate-basalganglia";}
            if ($filet eq "Brain-FrontalCortex(BA9)") {$filet = "Brain-FrontalCortex-BA9";}
            if ($filet eq "Brain-Nucleusaccumbens(basalganglia)") {$filet = "Brain-Nucleusaccumbens-basalganglia";}
            if ($filet eq "Brain-Putamen(basalganglia)") {$filet = "Brain-Putamen-basalganglia";}
            if ($filet eq "Skin-NotSunExposed(Suprapubic)") {$filet = "Skin-NotSunExposed-Suprapubic";}
            if ($filet eq "Skin-SunExposed(Lowerleg)") {$filet = "Skin-SunExposed-Lowerleg";}
            open (BETA, ">/group/im-lab/nas40t2/kaanan/PrediXcan/betas/allBetas/elasticNet${a}/TW_${filet}_exp_10-foldCV_elasticNet_hapmapSnpsCEU_all_chr1-22_unscaled.txt") or die "cant open allbetafile\n";
            print BETA "gene\trsid\tref\talt\tbeta\talpha\n";
            open (CV, ">/group/im-lab/nas40t2/kaanan/PrediXcan/betas/allResults/elasticNet${a}/TW_${filet}_unscaled.allResults.txt") or die "cant open cvfile\n";
            print CV "gene\talpha\tcvm\tlambda.iteration\tlambda.min\tn.snps\tR2\tpval\tgenename\n";
            foreach my $a (@alphas) {
                open (B, ">/group/im-lab/nas40t2/kaanan/PrediXcan/betas/allBetas/PrediXcanFormat/elasticNet${a}/TW_${filet}_exp_10-foldCV_elasticNet_alpha${a}_hapmapSnpsCEU_all_chr1-22_unscaled.txt") or die "cant open allbetafile\n";
                print B "gene\tSNP\trefAllele\teffectAllele\tbeta\n";
                ## example beta output
                # gene	rsid	ref	alt	beta	alpha
                # TUBB8	rs17159658	T	C	-1.072323e-02	0.05
                # TUBB8	rs11253269	T	C	-1.448800e-02	0.05
                # TUBB8	rs12761735	T	C	-6.895363e-03	0.05
                ## example results output
                # gene	alpha	cvm	lambda.iteration	lambda.min	n.snps	R2	pval
                # TUBB8	0.05	0.999109788111518	14	1.91159957460854	23	0.000900913386190969	0.581824483123699
                # ZMYND11	0.05	0.993386406424006	18	2.00048406597309	43	0.00536617027505977	0.17843845379816
                foreach my $chrom (@chromosomes) {
                    ## example beta input
                    # gene	SNP	refAllele	effectAllele	beta
                    # OR11H1	rs4819849	A	G	 0.0179785324
                    # OR11H1	rs1807512	T	C	-0.0272151097
                    # OR11H1	rs5994031	T	C	 0.0033450116
                    my $file = "/group/im-lab/nas40t2/kaanan/PrediXcan/betas/Resultsbychromosome/TW_${filet}_elasticNet_alpha${a}_hapmapSnpsCEU_weights_chr${chrom}_2015-11-12.txt";
                    open (IN, "$file") or die "cant open $file\n";
                    while (my $line = <IN>) {
                        chomp($line);
                        if ($line =~ /^gene/) {next;}
                        my @tmp = split(/\t/,$line);
                        print BETA "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t${a}\n";
                        print B "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\n";
                    }
                    close(IN);
                    ## example input file CV
                    # gene	alpha	cvm	lambda.iteration	lambda.min	n.snps	R2	pval
                    # POTEH	NA	NA	NA	NA	0	NA	NA
                    # OR11H1	0.5	0.999392341342453	6	0.212097780821214	12	0.00117975315933346	0.528527538066779
                    # CCT8L2	NA	NA	NA	NA	0	NA	NA
                    my $file2 = "/group/im-lab/nas40t2/kaanan/PrediXcan/betas/Resultsbychromosome/TW_${filet}_exp_10-foldCV_elasticNet_alpha${a}_hapmapSnpsCEU_chr${chrom}_2015-11-12.txt";
                    open (IN, "$file2") or die "cant open $file2\n";
                    while (my $line = <IN>) {
                        chomp($line);
                        if ($line =~ /^gene/) {next;}
                        print CV "$line\n";
                    }
                    close(IN);
                }
                close(B);
            }
            close(BETA);
            close(CV);
        }
    }
}
