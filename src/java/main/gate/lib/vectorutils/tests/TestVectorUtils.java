/*
 * Copyright (c) 2015-2016 The University Of Sheffield.
 *
 * This file is part of yodie-plugin-disambiguation. 
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
  
  VectorUtils vu;
  
  @Before
  public void prepare() {
    vu = new VectorUtils();
  }
  
  @Test
  public void testAdd() {
    double[] a = new double[]{1.0,2.0,3.0};
    double[] b = new double[]{2.0,1.1,3.0};
    double[] c = vu.add(a, b);
    assertEquals(c.length, a.length);
    assertEquals(3.0,c[0],vu.EPS);
    assertEquals(3.1,c[1],vu.EPS);
    assertEquals(6.0,c[2],vu.EPS);
  }
  
}
