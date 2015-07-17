/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.humangenome.util;

/**
 *
 * @author Yongzhuang Liu
 */
public class BaseUtils {

    public static char indexToBase(int index) {
        switch (index) {
            case 0:
                return 'A';
            case 1:
                return 'C';
            case 2:
                return 'G';
            case 3:
                return 'T';
            case 4:
                return 'N';
            default:
                return 0;
        }
    }

    public static int baseToIndex(char base) {
        base = Character.toUpperCase(base);
        switch (base) {
            case 'A':
                return 0;
            case 'C':
                return 1;
            case 'G':
                return 2;
            case 'T':
                return 3;
            case 'N':
                return 4;
            default:
                return -1;
        }
    }
}
