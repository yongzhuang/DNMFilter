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
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

/**
 *
 * @author Yongzhuang Liu
 */
public class PEDReader {

    String filename;

    public PEDReader(String filename){
        this.filename=filename; 
    }

    public List<Individual> getIndividuals() throws IOException {
        BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(filename)));
        String line = null;
        List<Individual> individuals = new ArrayList();
        while ((line = bufferedReader.readLine()) != null && line.trim().length() > 0) {
            StringTokenizer tokenizer = new StringTokenizer(line, "\n\t\" \"");
            String familyID = tokenizer.nextToken().trim();
            String individualID = tokenizer.nextToken().trim();
            String paternalID = tokenizer.nextToken().trim();
            String maternalID = tokenizer.nextToken().trim();
            int sex = Integer.parseInt(tokenizer.nextToken().trim());
            int phenotype = Integer.parseInt(tokenizer.nextToken().trim());
            individuals.add(new Individual(familyID, individualID, paternalID, maternalID, sex, phenotype));
        }
        bufferedReader.close();
        return individuals;
    }

    public List<Trio> getTrios() throws IOException {
        List<Individual> individuals = getIndividuals();
        List<Trio> trios = new ArrayList();
        for (int k = 0; k < individuals.size(); k = k + 3) {
            int paternalIndex = -1;
            int maternalIndex = -1;
            int offspringIndex = -1;
            for (int i = k; i < k + 3; i++) {
                Individual tmp = individuals.get(i);
                if (tmp.getPaternalID().equals("0") && tmp.getMaternalID().equals("0") && tmp.getSex() == 1) {
                    paternalIndex = i;
                } else if (tmp.getPaternalID().equals("0") && tmp.getMaternalID().equals("0") && tmp.getSex() == 2) {
                    maternalIndex = i;
                } else if (!tmp.getPaternalID().equals("0") && !tmp.getMaternalID().equals("0")) {
                    offspringIndex = i;
                } else {
                    break;
                }
            }
            if (paternalIndex != -1 && maternalIndex != -1 && offspringIndex != -1 && individuals.get(k).getFamilyID().equals(individuals.get(k + 1).getFamilyID())
                    && individuals.get(k + 1).getFamilyID().equals(individuals.get(k + 2).getFamilyID())) {
                trios.add(new Trio(individuals.get(k).getFamilyID(), individuals.get(paternalIndex), individuals.get(maternalIndex), individuals.get(offspringIndex)));
            }
        }
        return trios;
    }
    
}
