/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.humangenome.cmdline;

import edu.duke.humangenome.core.Extract;
import edu.duke.humangenome.core.GBM;
import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.Properties;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.log4j.Logger;

/**
 *
 * @author Yongzhuang Liu
 */
public class DNMFilter {

    private static Logger logger = Logger.getLogger(DNMFilter.class);

    public static void main(String[] args) throws IOException {
        String usage = "\nDNMFilter-0.1.0\n";
        usage = usage + "\nUsage: java -jar DNMFilter.jar <COMMAND> [OPTIONS]\n\n";
        usage = usage + "COMMANDS:\n"
                + "\textract\t\tExtract sequence features to build the training set\n"
                + "\tgbm\t\tUse gradient boosting approach to filter de novo mutations\n";
        String cmd = null;
        if (args.length > 0) {
            if (args[0].equals("extract") || args[0].equals("gbm")) {
                cmd = args[0];
            } else {
                logger.error("Command is not recognized!\n" + usage);
                return;
            }
        } else {
            System.out.println(usage);
            return;
        }
        run(cmd, args);
    }

    private static void run(String cmd, String[] args) throws IOException {
        long start = System.currentTimeMillis();
        CommandLineParser parser = new PosixParser();
        CommandLine commandLine = null;
        Options options = createOptions(cmd);
        try {
            if (options != null) {
                commandLine = parser.parse(options, args);
            }
        } catch (ParseException parseException) {
            logger.error("Invalid command line parameters!");
        }
        if (cmd.equals("extract")) {
            if (isValidated(commandLine, "extract")) {
                (new Extract(getProperties(commandLine, "extract"))).run();
            } else {
                printHelp(options,"extract");
                return;
            }
        }
        if (cmd.equals("gbm")) {
            if (isValidated(commandLine, "gbm")) {
                
               if(!(new GBM(getProperties(commandLine, "gbm"))).run()){
                   return;
               }
            } else {
                printHelp(options,"gbm");
                return;
            }
        }
        long end = System.currentTimeMillis();
        logger.info("Total running time is " + (end - start) / 1000 + " seconds");
        logger.info("Done!");
    }

    private static Options createOptions(String cmd) {
        Options options = new Options();
        if (cmd.equals("extract")) {
            options.addOption(OptionBuilder.withLongOpt("reference").withDescription("reference genome file (required)").hasArg().withArgName("FILE").create());
            options.addOption(OptionBuilder.withLongOpt("pedigree").withDescription("pedigree file (required)").hasArg().withArgName("FILE").create());
            options.addOption(OptionBuilder.withLongOpt("bam").withDescription("bam list file (required)").hasArg().withArgName("FILE").create());
            options.addOption(OptionBuilder.withLongOpt("output").withDescription("output file (required)").hasArg().withArgName("FILE").create());
            options.addOption(OptionBuilder.withLongOpt("positive").withDescription("known true positive DNM file (required)").hasArg().withArgName("FILE").create());
            options.addOption(OptionBuilder.withLongOpt("negative").withDescription("known false positive DNM file (required)").hasArg().withArgName("FILE").create());
            return options;
        } else if (cmd.equals("gbm")) {
            options.addOption(OptionBuilder.withLongOpt("reference").withDescription("reference genome file (required)").hasArg().withArgName("FILE").create());
            options.addOption(OptionBuilder.withLongOpt("pedigree").withDescription("pedigree file (required)").hasArg().withArgName("FILE").create());
            options.addOption(OptionBuilder.withLongOpt("bam").withDescription("bam list file (required)").hasArg().withArgName("FILE").create());
            options.addOption(OptionBuilder.withLongOpt("output").withDescription("output file (required)").hasArg().withArgName("FILE").create());
            options.addOption(OptionBuilder.withLongOpt("candidate").withDescription("candidate DNM file (required)").hasArg().withArgName("FILE").create());
            options.addOption(OptionBuilder.withLongOpt("training").withDescription("training set (required)").hasArg().withArgName("FILE").create());
            options.addOption(OptionBuilder.withLongOpt("configuration").withDescription("feature configuration file (required)").hasArg().withArgName("FILE").create());
            options.addOption(OptionBuilder.withLongOpt("cutoff").withDescription("cutoff to determine a putative DNM (optional, default 0.4)").hasArg().withArgName("DOUBLE").create());
            return options;
        } else {
            return null;
        }
    }

    private static Properties getProperties(CommandLine line, String cmd) {
        Properties properties = new Properties();
        if (cmd.equals("extract")) {
            properties.put("reference", line.getOptionValue("reference"));
            properties.put("pedigree", line.getOptionValue("pedigree"));
            properties.put("bam", line.getOptionValue("bam"));
            properties.put("output", line.getOptionValue("output"));
            properties.put("positive", line.getOptionValue("positive"));
            properties.put("negative", line.getOptionValue("negative"));
        }
        if (cmd.equals("gbm")) {
            properties.put("reference", line.getOptionValue("reference"));
            properties.put("pedigree", line.getOptionValue("pedigree"));
            properties.put("bam", line.getOptionValue("bam"));
            properties.put("output", line.getOptionValue("output"));
            properties.put("candidate", line.getOptionValue("candidate"));
            properties.put("training", line.getOptionValue("training"));
            properties.put("configuration", line.getOptionValue("configuration"));
            if (!line.hasOption("cutoff")) {
                properties.put("cutoff", "0.4");    
            }
            else{
                properties.put("cutoff", line.getOptionValue("cutoff"));
            }
        }
        return properties;
    }

    private static boolean isValidated(CommandLine line, String cmd) {
        boolean tag = true;
        if (!line.hasOption("reference") || !(new File(line.getOptionValue("reference")).isFile())) {
            logger.error("The reference genome file is not correctly specified!");
            tag = false;
        }
        if (!line.hasOption("pedigree") || !(new File(line.getOptionValue("pedigree")).isFile())) {
            logger.error("The pedigree file is not correctly specified!");
            tag = false;
        }
        if (!line.hasOption("bam") || !(new File(line.getOptionValue("bam")).isFile())) {
            logger.error("The bam list file is not correctly specified!");
            tag = false;
        }
        if (cmd.equals("extract")) {
            if (!line.hasOption("positive") || !(new File(line.getOptionValue("positive")).isFile())) {
                logger.error("The known true positive DNM file is not correctly specified!");
                tag = false;
            }
            if (!line.hasOption("negative") || !(new File(line.getOptionValue("negative")).isFile())) {
                logger.error("The known false positive DNM file is not correctly specified!");
                tag = false;
            }
        }
        if (cmd.equals("gbm")) {
            if (!line.hasOption("training") || !(new File(line.getOptionValue("training")).isFile())) {
                logger.error("The training set is not correctly specified!");
                tag = false;
            }
            if (!line.hasOption("candidate") || !(new File(line.getOptionValue("candidate")).isFile())) {
                logger.error("The candidate DNM file is not correctly specified!");
                tag = false;
            }
            if (!line.hasOption("configuration") || !(new File(line.getOptionValue("configuration")).isFile())) {
                logger.error("The feature configuration file is not correctly specified!");
                tag = false;
            }
        }
        if (!line.hasOption("output")) {
            logger.error("The output file is not correctly specified!");
            tag = false;
        }
        return tag;
    }

    private static void printHelp(Options options,String command) {
        System.out.println();
        String cmdLineSyntax = "java -jar DNMFilter.jar "+command+" [OPTIONS]\n";
        HelpFormatter formatter = new HelpFormatter();
        formatter.setOptionComparator(new OptionComarator());
        formatter.printHelp(cmdLineSyntax, options);
    }

    private static class OptionComarator<T extends Option> implements Comparator<T> {

        private static final String OPTS_ORDER = "rpnbktco";

        public int compare(T o1, T o2) {
            return OPTS_ORDER.indexOf(o1.getLongOpt().charAt(0)) - OPTS_ORDER.indexOf(o2.getLongOpt().charAt(0));
        }
    }
}
