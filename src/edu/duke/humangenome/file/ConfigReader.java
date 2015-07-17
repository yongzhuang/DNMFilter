/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.humangenome.file;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Properties;
import java.util.StringTokenizer;

/**
 *
 * @author Yongzhuang Liu
 */
public class ConfigReader {

    private String filename;

    public ConfigReader(String filename){
        this.filename=filename;
    }

    public Properties parse() throws IOException {
        BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(filename)));
        Properties properties = new Properties();
        String line = null;
        while ((line = bufferedReader.readLine()) != null && line.trim().length() > 0) {
            StringTokenizer tokenizer = new StringTokenizer(line, "\n\t\" \"");
            String sampleID = tokenizer.nextToken().trim();
            String sampleFile = tokenizer.nextToken().trim();
            properties.put(sampleID, sampleFile);
        }
        bufferedReader.close();
        return properties;
    }

}
