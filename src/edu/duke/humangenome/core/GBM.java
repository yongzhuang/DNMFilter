/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.humangenome.core;

import edu.duke.humangenome.file.ConfigReader;
import edu.duke.humangenome.file.DNMReader;
import edu.duke.humangenome.file.DNMRecord;
import edu.duke.humangenome.file.PEDReader;
import edu.duke.humangenome.file.Trio;
import edu.duke.humangenome.sam.MultiPileup;
import edu.duke.humangenome.sam.MultiSamLocusIterator;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFileWalker;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.picard.util.SamLocusIterator;
import net.sf.samtools.SAMFileReader;
import org.apache.log4j.Logger;
import rcaller.RCaller;
import rcaller.RCode;

/**
 *
 * @author Yongzhuang Liu
 */
public class GBM {

    private static Logger logger = Logger.getLogger(GBM.class);
    private final static int FATHER_INDEX = 0;
    private final static int MOTHER_INDEX = 1;
    private final static int OFFSPRING_INDEX = 2;
    private Properties properties;

    public GBM(Properties properties) {
        this.properties = properties;
    }

    public boolean run() throws IOException {
        
        logger.info("Checking R environment and packages ......");
 
        Process pid = Runtime.getRuntime().exec("which Rscript");
        BufferedReader runTimeReader = new BufferedReader(new InputStreamReader(pid.getInputStream()));
        String R_HOME = runTimeReader.readLine().trim();
        if (R_HOME.equals("")) {
            logger.error("Rscript exectuable is not set in the PATH environment variable!");
            return false;
        }
        
        if (!checkInstalledPackages(R_HOME)) {
            return false;
        }

        logger.info("Extracting features ......");
        
        int numOfAllSites=0;
        String outputPath = (new File(properties.getProperty("output"))).getAbsolutePath();
        String folder = outputPath.substring(0, outputPath.lastIndexOf("/"));
        File tmpTesting = File.createTempFile("Testing", ".csv", new File(folder));
        tmpTesting.deleteOnExit();
        
        PrintWriter testing = new PrintWriter(tmpTesting);
        testing.write("Family_ID" + "," + "Chrom" + "," + "Position" + "," );
        testing.write("Father_Allele_Balance" + "," + "Father_Mean_Base_Quality_For_Ref" + "," + "Father_Mean_Base_Quality_For_Alt" + "," + "Father_Read_Depth" + ",");
        testing.write("Father_Mean_Mapping_Quality_For_Ref" + "," + "Father_Mean_Mapping_Quality_For_Alt" + "," + "Father_Mean_Distance_To_Three_Prime_For_Ref" + "," + "Father_Mean_Distance_To_Three_Prime_For_Alt" + ",");
        testing.write("Father_Fraction_Of_MQ0_Reads_For_Ref" + "," + "Father_Fraction_Of_MQ0_Reads_For_Alt" + ",");
        testing.write("Father_Fraction_Of_Soft_Clipped_Reads_For_Ref" + "," + "Father_Fraction_Of_Soft_Clipped_Reads_For_Alt" + ",");
        testing.write("Father_Mean_Nearby_Mismatches_For_Ref" + "," + "Father_Mean_Nearby_Mismatches_For_Alt" + ",");
        testing.write("Father_Mean_Nearby_Indels_For_Ref" + "," + "Father_Mean_Nearby_Indels_For_Alt" + ",");

        testing.write("Father_Strand_Direction_For_Ref" + "," + "Father_Strand_Direction_For_Alt" + ",");
        testing.write("Father_Strand_Bias" + ",");

        testing.write("Mother_Allele_Balance" + "," + "Mother_Mean_Base_Quality_For_Ref" + "," + "Mother_Mean_Base_Quality_For_Alt" + "," + "Mother_Read_Depth" + ",");
        testing.write("Mother_Mean_Mapping_Quality_For_Ref" + "," + "Mother_Mean_Mapping_Quality_For_Alt" + "," + "Mother_Mean_Distance_To_Three_Prime_For_Ref" + "," + "Mother_Mean_Distance_To_Three_Prime_For_Alt" + ",");
        testing.write("Mother_Fraction_Of_MQ0_Reads_For_Ref" + "," + "Mother_Fraction_Of_MQ0_Reads_For_Alt" + ",");
        testing.write("Mother_Fraction_Of_Soft_Clipped_Reads_For_Ref" + "," + "Mother_Fraction_Of_Soft_Clipped_Reads_For_Alt" + ",");
        testing.write("Mother_Mean_Nearby_Mismatches_For_Ref" + "," + "Mother_Mean_Nearby_Mismatches_For_Alt" + ",");
        testing.write("Mother_Mean_Nearby_Indels_For_Ref" + "," + "Mother_Mean_Nearby_Indels_For_Alt" + ",");
        testing.write("Mother_Strand_Direction_For_Ref" + "," + "Mother_Strand_Direction_For_Alt" + ",");
        testing.write("Mother_Strand_Bias" + ",");

        testing.write("Offspring_Allele_Balance" + "," + "Offspring_Mean_Base_Quality_For_Ref" + "," + "Offspring_Mean_Base_Quality_For_Alt" + "," + "Offspring_Read_Depth" + ",");
        testing.write("Offspring_Mean_Mapping_Quality_For_Ref" + "," + "Offspring_Mean_Mapping_Quality_For_Alt" + "," + "Offspring_Mean_Distance_To_Three_Prime_For_Ref" + "," + "Offspring_Mean_Distance_To_Three_Prime_For_Alt" + ",");
        
        testing.write("Offspring_Fraction_Of_MQ0_Reads_For_Ref" + "," + "Offspring_Fraction_Of_MQ0_Reads_For_Alt" + ",");
        testing.write("Offspring_Fraction_Of_Soft_Clipped_Reads_For_Ref" + "," + "Offspring_Fraction_Of_Soft_Clipped_Reads_For_Alt" + ",");
        testing.write("Offspring_Mean_Nearby_Mismatches_For_Ref" + "," + "Offspring_Mean_Nearby_Mismatches_For_Alt" + ",");
        testing.write("Offspring_Mean_Nearby_Indels_For_Ref" + "," + "Offspring_Mean_Nearby_Indels_For_Alt" + ",");
        testing.write("Offspring_Strand_Direction_For_Ref" + "," + "Offspring_Strand_Direction_For_Alt" + ",");
        testing.write("Offspring_Strand_Bias" + ",");

        testing.write("PValue_Father_To_Offspring" + "," + "PValue_Mother_To_Offspring" + "\n");

        File referenceSequenceFile = new File(properties.getProperty("reference"));
        List<Trio> trios = (new PEDReader(properties.getProperty("pedigree"))).getTrios();
        Properties bams = (new ConfigReader(properties.getProperty("bam"))).parse();
        PrintWriter output = new PrintWriter(properties.getProperty("output"));
        DNMReader dnmReader = new DNMReader(properties.getProperty("candidate"));
        double cutoff=Double.parseDouble(properties.getProperty("cutoff"));
        HashMap<String, Trio> map = new HashMap();
        for (int i = 0; i < trios.size(); i++) {
            map.put(trios.get(i).getFamilyID(), trios.get(i));
        }
        
        DNMReader countReader = new DNMReader(properties.getProperty("candidate"));
        DNMRecord dnmRecord=null;
        int count=0;
        while ((dnmRecord=countReader.getNextRecord()) != null) {
            count++;
        }
        countReader.close();
        
        dnmRecord = dnmReader.getNextRecord();
        String familyID = null;
        int numOfTrios=0;
        while (dnmRecord != null) {
            numOfTrios++;
            familyID = dnmRecord.getFamilyID();
            Trio trio = map.get(familyID);
            MultiSamLocusIterator trioSamLocusIterator = new MultiSamLocusIterator();
            SAMFileReader[] trioSAMFileReader = new SAMFileReader[3];
            IntervalList[] trioIntervalList = new IntervalList[3];
            SamLocusIterator[] samLocusIterator = new SamLocusIterator[3];
            trioSAMFileReader[0] = new SAMFileReader(new File(bams.getProperty(trio.getFather().getIndividualID())));
            trioSAMFileReader[1] = new SAMFileReader(new File(bams.getProperty(trio.getMother().getIndividualID())));
            trioSAMFileReader[2] = new SAMFileReader(new File(bams.getProperty(trio.getOffspring().getIndividualID())));
            for (int i = 0; i < 3; i++) {
                trioSAMFileReader[i].setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
                trioIntervalList[i] = new IntervalList(trioSAMFileReader[i].getFileHeader());
            }
            while (dnmRecord != null && dnmRecord.getFamilyID().equals(familyID)) {
                Interval interval = new Interval(dnmRecord.getChrom(), dnmRecord.getPos(), dnmRecord.getPos());
                for (int i = 0; i < 3; i++) {
                    trioIntervalList[i].add(interval);
                }
                dnmRecord = dnmReader.getNextRecord();
            }
            for (int i = 0; i < 3; i++) {
                //trioIntervalList[i].sort();
                samLocusIterator[i] = new SamLocusIterator(trioSAMFileReader[i], trioIntervalList[i], true);
            }
            trioSamLocusIterator.add(samLocusIterator[0], FATHER_INDEX);
            trioSamLocusIterator.add(samLocusIterator[1], MOTHER_INDEX);
            trioSamLocusIterator.add(samLocusIterator[2], OFFSPRING_INDEX);
            ReferenceSequence referenceSequence = null;
            Map<Integer, SamLocusIterator.LocusInfo> tmp = null;
            ReferenceSequenceFileWalker referenceSequenceFileWalker = new ReferenceSequenceFileWalker(referenceSequenceFile);
            while ((tmp = trioSamLocusIterator.getLocusInfos()) != null) {
                MultiPileup trioPileup = new MultiPileup(tmp);
                if (referenceSequence == null || referenceSequence.getContigIndex() != trioPileup.getReferenceIndex()) {
                    referenceSequence = referenceSequenceFileWalker.get(trioPileup.getReferenceIndex());
                }
                FeatureSelection featureSelection = new FeatureSelection(referenceSequence, trioPileup);
                testing.write(familyID + "," + trioPileup.getReferenceName() + "," + trioPileup.getPosition() + "," + featureSelection.extract() + "\n");
                numOfAllSites++;
            }
            for (int i = 0; i < 3; i++) {
                trioSAMFileReader[i].close();
            }
            trioSamLocusIterator.close();
        }
        dnmReader.close();
        testing.close();
        
        logger.info("Filtering DNMs ......");
        
        if (numOfAllSites == 0) {
            logger.info((count - numOfAllSites) + " non-informative candidate DNMs were skipped.");
            logger.info(numOfAllSites + " informative candidate DNMs in " + numOfTrios + " trios were processed.");
            logger.info("Number of true DNMs is " + 0 + "," + " Number of false positive DNMs is " + 0 + ".");
            return true;
        }
    
        Properties features = new Properties();
        features.load(new FileReader(properties.getProperty("configuration")));
        List<String> selectedFeaturesList=new ArrayList();
        
        for(Object key : features.keySet()){
            int value=Integer.parseInt(((String)features.get(key)).trim());
            if(value==1){
                selectedFeaturesList.add((String)key);
            }  
        }
        String[] selectedFeatures=new String[selectedFeaturesList.size()];
        selectedFeaturesList.toArray(selectedFeatures);
        
        String trainingPath = properties.getProperty("training");
        double[] result=null;
        try {
            
            RCaller caller = new RCaller();
            RCode code = new RCode();
            caller.setRscriptExecutable(R_HOME);
            String[] input = new String[2];
            input[0] = trainingPath;
            input[1] = tmpTesting.getAbsolutePath();
            code.addStringArray("input", input);
            code.addStringArray("selectedFeatures", selectedFeatures);
            code.addRCode("tryCatch({");
            code.addRCode("require(Runiversal)");
            code.addRCode("require(gbm)");
            code.addRCode("training<-read.table(input[1],header=T, sep = \",\")");
            code.addRCode("testing<-read.table(input[2],header=T, sep = \",\")");
            code.addRCode("training$Class<- ifelse(training$Class==\"positive\",1,0)");
            code.addRCode("gbm <- gbm(Class~.,data=training[,c(selectedFeatures,\"Class\")],shrinkage=0.001,distribution='bernoulli',cv.folds=10,n.trees=50000,interaction.depth = 1, bag.fraction=0.5, verbose=F)");
            code.addRCode("best.iter <- gbm.perf(gbm,method='cv',plot.it=FALSE)");
            code.addRCode("pred<-predict(gbm, testing[,selectedFeatures], best.iter, type=\"response\" )");
            code.addRCode("}, error=function(err){");
            code.addRCode("print(err)");
            code.addRCode("})");
            caller.setRCode(code);
            caller.runAndReturnResult("pred");
            result = caller.getParser().getAsDoubleArray("pred");
        } catch (Exception e) {
        
        }
        
        BufferedReader bufferedReader = new BufferedReader(new FileReader(tmpTesting));
        String line;
        int k = 0;
        line = bufferedReader.readLine();
        int numOfTrueDNMs = 0;
        while ((line = bufferedReader.readLine()) != null) {
            String[] record = line.split(",");
            if (result[k] >= cutoff) {
                numOfTrueDNMs++;
                output.write(record[0] + "," + record[1] + "," + record[2] + "," + result[k] + "\n");
            }
            k++;
        }
        output.close();
        logger.info((count-numOfAllSites)+" non-informative candidate DNMs were skipped.");
        logger.info(numOfAllSites+" informative candidate DNMs in " + numOfTrios+" trios were processed.");
        logger.info("Number of true DNMs is "+numOfTrueDNMs+","+" Number of false positive DNMs is "+(numOfAllSites-numOfTrueDNMs)+".");
        return true;
    }
    
    
    public boolean checkInstalledPackages(String R_HOME){
        
            boolean tag=true;
            RCaller caller = new RCaller();
            caller.setRscriptExecutable(R_HOME);
            caller.cleanRCode();
            RCode code = new RCode();
            String[] packages=new String[]{"Runiversal","gbm"};
            code.addStringArray("packages", packages);
            code.addRCode("label<-c(1, 1)");
            code.addRCode("for(i in 1:2){");
            code.addRCode("if(!require(packages[i], character.only=TRUE)){");
            code.addRCode("label[i]=0");
            code.addRCode("}");
            code.addRCode("}");
            caller.setRCode(code);
            caller.runAndReturnResult("label");
            int[] label=caller.getParser().getAsIntArray("label");
            for(int i=0;i<packages.length;i++){
                if(label[i]==0){
                logger.error(packages[i]+" is not installed in R environment!" );
                tag=false; 
                }
            }
            return tag;
    }
}
