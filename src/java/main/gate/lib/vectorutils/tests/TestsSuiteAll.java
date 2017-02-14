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

import org.junit.runner.RunWith;
import org.junit.runners.Suite;

/**
 * Suite for all the tests.
 * @author Johann Petrak
 */

@RunWith(Suite.class)
@Suite.SuiteClasses({
  TestVectorUtils.class,
})
public class TestsSuiteAll {
  public static void main(String args[]) {
    org.junit.runner.JUnitCore.main(TestsSuiteAll.class.getCanonicalName());
  }  
  
}
