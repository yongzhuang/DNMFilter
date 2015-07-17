/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.humangenome.file;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 *
 * @author Yongzhuang Liu
 */
public class DNMReader {

    private BufferedReader bufferedReader;

    public DNMReader(String filename) throws IOException {
        this.bufferedReader = new BufferedReader(new FileReader(new File(filename)));
    }

    public DNMRecord getNextRecord() throws IOException {
        String line;
        while ((line = bufferedReader.readLine()) != null && line.trim().length() > 0) {
            String[] record = line.split(",");
            String familyID = record[0];
            String chrom = record[1];
            int pos = Integer.parseInt(record[2]);
            if (record.length == 3) {
                return new DNMRecord(familyID, chrom, pos, null, null, null);
            } else if (record.length == 4) {
                String ref = record[3];
                return new DNMRecord(familyID, chrom, pos, ref, null, null);
            } else if (record.length == 4) {
                String ref = record[3];
                String var = record[4];
                return new DNMRecord(familyID, chrom, pos, ref, var, null);
            } else {
                String ref = record[3];
                String var = record[4];
                String[] restOfFields = new String[record.length - 5];
                for (int i = 5; i < record.length; i++) {
                    restOfFields[i - 5] = record[i];
                }
                return new DNMRecord(familyID, chrom, pos, ref, var, restOfFields);
            }
        }
        return null;
    }

    public void close() throws IOException {
        bufferedReader.close();
    }
}
