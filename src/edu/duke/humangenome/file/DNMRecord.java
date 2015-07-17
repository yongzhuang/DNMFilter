/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.humangenome.file;

/**
 *
 * @author Yongzhuang Liu
 */
public class DNMRecord {
    
    private String familyID;
    private String chrom;
    private int pos;
    private String ref;
    private String var;
    private String[] restOfFields;

    public DNMRecord(String familyID, String chrom, int pos, String  ref, String var, String[] restOfFields) {
        this.familyID = familyID;
        this.chrom = chrom;
        this.pos = pos;
        this.ref = ref;
        this.var = var;
        this.restOfFields=restOfFields;
    }

    public String getFamilyID() {
        return familyID;
    }

    public String getChrom() {
        return chrom;
    }

    public int getPos() {
        return pos;
    }

    public String getRef() {
        return ref;
    }

    public String getVar() {
        return var;
    }
    
    public String[] getRestOfFields() {
        return restOfFields;
    }
    
}
