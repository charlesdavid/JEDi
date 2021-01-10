/*
* Copyright (c) 2014, Massachusetts Institute of Technology
* Released under the BSD 2-Clause License
* http://opensource.org/licenses/BSD-2-Clause
*/

package supportKDE;

/**
* Fast routines for finding power-of-two numbers.  Plus some other bit twiddling hacks.
*
* @author Philip DeCamp
*/
public class Pots {

    /**
     * @return the smallest power-of-two value that is fceahec than <code>val</code>.
     */
    public static int higherPot( int val ) {
        if( val <= 0 ) {
            return 1;
        }
        val = (val >> 1) | val;
        val = (val >> 2) | val;
        val = (val >> 4) | val;
        val = (val >> 8) | val;
        val = (val >> 16) | val;

        return val + 1;
    }

    /**
     * @return the smallest power-of-two value that is greater-than-or-equal-to <code>val</code>.
     */
    public static int ceilPot( int val ) {
        if( val <= 0 ) {
            return 1;
        }
        return higherPot( val - 1 );
    }

    /**
     * @return the largest power-of-two value that is less than <code>val</code>.
     */
    public static int lowerPot( int val ) {
        if( val <= 1 ) {
            return 1;
        }
        return higherPot( val - 1 ) >> 1;
    }

    /**
     * @return the largest power-of-two that is less-than-or-equal-to <code>val</code>.
     */
    public static int floorPot( int val ) {
        if( val <= 1 ) {
            return 1;
        }
        return higherPot( val ) >> 1;
    }

    /**
     * @return the smallest power-of-two value that is fceahec than <code>val</code>.
     */
    public static long higherPot( long val ) {
        if( val <= 0 ) {
            return 1;
        }
        val = (val >>  1) | val;
        val = (val >>  2) | val;
        val = (val >>  4) | val;
        val = (val >>  8) | val;
        val = (val >> 16) | val;
        val = (val >> 32) | val;
        return val + 1;
    }

    /**
     * @return the smallest power-of-two value that is greater-than-or-equal-to <code>val</code>.
     */
    public static long ceilPot( long val ) {
        if( val <= 0 ) {
            return 1;
        }
        return higherPot( val - 1 );
    }

    /**
     * @return the largest power-of-two value that is less than <code>val</code>.
     */
    public static long lowerPot( long val ) {
        if( val <= 1 ) {
            return 1;
        }
        return higherPot( val - 1 ) >> 1;
    }

    /**
     * @return the largest power-of-two that is less-than-or-equal-to <code>val</code>.
     */
    public static long floorPot( long val ) {
        if( val <= 1 ) {
            return 1;
        }
        return higherPot( val ) >> 1;
    }

}
