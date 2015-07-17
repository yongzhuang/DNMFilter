/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.humangenome.file;

/**
 *
 * @author Yongzhuang Liu
 */
public class Individual {

    private String familyID;
    private String individualID;
    private String paternalID;
    private String maternalID;
    private int sex;
    private int phenotype;

    public Individual(String familyID, String individualID, String paternalID, String maternalID, int sex, int phenotype) {
        this.familyID = familyID;
        this.individualID = individualID;
        this.paternalID = paternalID;
        this.maternalID = maternalID;
        this.sex = sex;
        this.phenotype = phenotype;
    }

    public String getFamilyID() {
        return familyID;
    }

    public String getIndividualID() {
        return individualID;
    }

    public String getMaternalID() {
        return maternalID;
    }

    public String getPaternalID() {
        return paternalID;
    }

    public int getSex() {
        return sex;
    }

    public int getPhenotype() {
        return phenotype;
    }
}
