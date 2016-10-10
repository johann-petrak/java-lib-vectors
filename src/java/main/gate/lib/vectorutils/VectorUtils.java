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
package gate.lib.vectorutils;

import java.util.Arrays;
import org.apache.log4j.Logger;
import java.util.List;

/**
 * Simple collection of vector utilities. 
 * 
 * @author Johann Petrak
 */
public class VectorUtils {
  
  public double EPS = 10.0*Double.MIN_VALUE;
  
  // If true, throws an exception if some error occurs, otherwise logs
  // an ERROR and continues with some fake value.
  boolean throwOnError = true;
  
  /**
   * Default constructor, creates the instance with the default option
   * settings.
   * 
   */
  public VectorUtils() {
  }
  
  /**
   * Create the instance and set the options.
   * @param eps
   * @param throwOnError 
   */
  public VectorUtils(double eps, boolean throwOnError) {
    EPS = eps;
    this.throwOnError = throwOnError;
  }
  
  Logger logger = Logger.getLogger(VectorUtils.class);
  
  /////// SINGLE VECTOR UTILS
  
  /**
   * Normalize the vector by its L2 norm.
   * @param vector
   * @return a copy of the vector, normalized by the L2 norm of the original.
   */
  public double[] normalizeL2(double[] vector) {
    double norm = normL2(vector);
    if(isZero(norm)) {
      if(throwOnError)
        throw new RuntimeException("Cannot normalize vector, L2 norm too close to zero");
      else {
        logger.error("normalizeL2: norm too close to zero, using 1.0");
        return Arrays.copyOf(vector, vector.length);        
      }
    } 
    double[] ret = new double[vector.length];
    for(int i=0; i<vector.length; i++) {
      ret[i] = vector[i] / norm;
    }
    return ret;
  }

  /**
   * Normalize the vector by its L1 norm.
   * @param vector
   * @return a copy of the vector, normalized by the L1 norm of the original.
   */
  public double[] normalizeL1(double[] vector) {
    double norm = normL1(vector);
    if(isZero(norm)) {
      throw new RuntimeException("Cannot normalize vector, L2 norm too close to zero");
    }
    double[] ret = new double[vector.length];
    for(int i=0; i<vector.length; i++) {
      ret[i] = vector[i] / norm;
    }
    return ret;
  }
  
  /////// functions on two or more vectors, returning a scalar
  
  /**
   * Calculate the L2 norm of the vector
   * @param vector
   * @return 
   */
  public double normL2(double[] vector) {
    double sumsquares = 0.0;
    for(double val : vector) {
      sumsquares += val * val;
    }
    double norm = Math.sqrt(sumsquares);
    return norm;
  }

  /**
   * Calculate the L1 norm of the vector
   * @param vector
   * @return 
   */
  public double normL1(double[] vector) {
    double norm = 0.0;
    for(double val : vector) {
      norm += Math.abs(val);
    }
    return norm;
  }
  
  /**
   * Calculate the dot-product (inner product) of two vectors.
   * 
   * @param vector1
   * @param vector2
   * @return dot-product of the vectors
   */
  public double dot(double[] vector1, double[] vector2) {
    if (vector1.length != vector2.length) {
      throw new RuntimeException("Cannot add vectors of different size");
    }
    double ret = 0.0;
    for(int i=0; i<vector1.length; i++) {
      ret += vector1[i] * vector2[i];
    }
    return ret;
  }
  
  /**
   * Calculate the cosine similarity of two vectors.
   * 
   * @param vector1
   * @param vector2
   * @return cosine-similarity of the two vectors
   */
  public double simCosine(double[] vector1, double[] vector2) {
    if(vector1.length != vector2.length) {
      throw new RuntimeException("Vectos must have equal length for simCosine");
    }
    double div = normL2(vector1) * normL2(vector2);
    // TODO: what if div is zero
    return dot(vector1,vector2) / div;
  }
  
  
  /////// functions on two vectors or lists of vectors which create a 
  /////// single vector result
  
  /**
   * Return a vector that is the sum of all the vectors passed in the list.
   * Each dimension in the result vector is the sum of the same dimension
   * from all vectors in the list.
   * 
   * @param vectors
   * @return sum vector
   */
  public double[] sum(List<double[]> vectors) {
    return sum(vectors, null);
  }
  
  
  /**
   * Return a vector that is the weighted sum of all the vectors passed in the list.
   * Each dimension in the result vector is the sum of the same dimension
   * from all vectors in the list, weighted by the weight for the corresponding vector.
   * 
   * If the weights are null, no weighting is done and a normal sum is calculated.
   * 
   * @param vectors
   * @return sum vector
   */
  public double[] sum(List<double[]> vectors, List<Double> weights) {
    if(vectors.isEmpty()) {
      throw new RuntimeException("Cannot create sum of no vectors");
    }
    if(weights != null && weights.size() != vectors.size()) {
      throw new RuntimeException("Weights must have same size as vectors");
    }
    int len = vectors.get(0).length;
    double[] ret = new double[len];
    for(int i=0;i<vectors.size();i++) {
      double[] vec = vectors.get(i);
      if(vec.length != len) {
        throw new RuntimeException("Vectors for sum of different lengths: "+len+"/"+vec.length);
      }
      for(int j=0;j<len;j++) {
        if(weights != null)
          ret[j] = ret[j] + (vec[j] * weights.get(i));
        else 
          ret[j] = ret[j] + vec[j];
      }
    }
    return ret;
  }
  
  /**
   * Return a vector that is the weighted average of all the vectors passed in the list.
   * Each dimension in the result vector is the average of the same dimension
   * from all vectors in the list, weighted by the weight for the corresponding vector.
   * 
   * If the weights are null, no weighting is done and an unweighted average is calculated.
   * 
   * @param vectors
   * @return sum vector
   */
  public double[] average(List<double[]> vectors, List<Double> weights) {
    double[] ret = sum(vectors,weights);
    double sum = 0;
    if(weights == null) {
      sum = vectors.size();
    } else {
      for(Double weight : weights) {
        sum += weight;
      }
      if(isZero(sum)) {
        throw new RuntimeException("Sum of weights is zero or too close to zero");
      }
    }
    for(int i=0; i<ret.length; i++) {
      ret[i] = ret[i] / sum;
    }
    return ret;
  }

  /**
   * Calculate a vector where each dimension is the maximum of that dimension
   * among all the given vectors. 
   * @param vectors
   * @return maxpooling vector
   */
  public double[] maxpooling(List<double[]> vectors) {
    if(vectors.isEmpty()) {
      throw new RuntimeException("Cannot create maximum pooling of no vectors");
    }
    int len = vectors.get(0).length;
    double[] ret = new double[len];
    for (double[] vec : vectors) {
      if(vec.length != len) {
        throw new RuntimeException("Vectors for sum of different lengths: "+len+"/"+vec.length);
      }
      // For max pooling, replace each dimension by the current vectors value
      // if that is larger. 
      for(int j=0;j<len;j++) {
        if(vec[j] > ret[j])
          ret[j] = vec[j];
      }
    }
    return ret;
  }
  
  /**
   * Return a vector that is the sum of the two given vectors
   * @param vec1
   * @param vec2
   * @return 
   */
  public double[] add(double[] vec1, double[] vec2) {
    if (vec1.length != vec2.length) {
      throw new RuntimeException("Cannot add vectors of different size");
    } else {
      double[] ret = new double[vec1.length];
      for (int i = 0; i < vec1.length; i++) {
        ret[i] = vec1[i] + vec2[2];
      }
      return ret;
    }
  }
  
  /**
   * Return a vector that is the element-wise product of the two vectors.
   * @param vector1
   * @param vector2
   * @return 
   */
  public double[] mul(double[] vector1, double[] vector2) {
    if (vector1.length != vector2.length) {
      throw new RuntimeException("Cannot add vectors of different size");
    }
    double[] ret = new double[vector1.length];
    for(int i=0; i<vector1.length; i++) {
      ret[i] = vector1[i] * vector2[i];
    }
    return ret;
  }
  
  public boolean isZero(double value) {
    return (Math.abs(value) < EPS);
  }
  
  
}
