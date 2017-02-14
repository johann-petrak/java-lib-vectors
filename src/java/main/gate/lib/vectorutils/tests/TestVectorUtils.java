/*
 * Copyright (c) 2015-2016 The University Of Sheffield.
 *
 * This file is part of java-lib-vectors
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 2.1 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software. If not, see <http://www.gnu.org/licenses/>.
 */
package gate.lib.vectorutils.tests;

import gate.lib.vectorutils.VectorUtils;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import org.junit.BeforeClass;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import org.junit.Before;

/**
 * Tests for the VectorUtils class.
 * @author Johann Petrak
 */
public class TestVectorUtils {

  @BeforeClass
  public static void init() {
    
  }
  
  double LEPS = 0.000001;
  
  VectorUtils vu;
  
  @Before
  public void prepare() {
    vu = new VectorUtils();
  }
  
  @Test
  public void test1() {
    double[] a = new double[]{1.0,2.0,3.0};
    double[] b = new double[]{2.0,1.1,3.0};
    double[] t_add = vu.add(a, b);
    assertEquals(t_add.length, a.length);
    assertEquals(3.0,t_add[0],vu.EPS);
    assertEquals(3.1,t_add[1],vu.EPS);
    assertEquals(6.0,t_add[2],vu.EPS);
    
    assertEquals(13.2,vu.dot(a, b),vu.EPS);
    
    assertEquals(3.7416573867739413,vu.normL2(a),LEPS);
    
    // we got the number from scipy.spatial.distance.cosine which is a distance,
    // so have to convert it back into a similarity
    assertEquals((1.0-0.06413571368455151),vu.simCosine(a, b),LEPS);
    
    double[] c = new double[]{2.2,0.3,1.1};
    
    assertEquals((1.0-0.3420670580023687),vu.simCosine(a, c), LEPS);
    
    assertEquals((1.0-1.5765566601970549),vu.simCorrelation(a, c), LEPS);
    
    assertEquals(2.817800560721074,vu.distEuclidean(a, c),LEPS);
    
    List<double[]> vectors = new ArrayList<double[]>();
    vectors.add(a);
    vectors.add(b);
    vectors.add(c);
    
    double avg[] = vu.average(vectors);
    assertEquals(1.73333333,avg[0],LEPS);
    assertEquals(1.13333333,avg[1],LEPS);
    assertEquals(2.36666667,avg[2],LEPS);
    
  }
  
}
