package com.csc205.project1;

import java.util.logging.Logger;
import java.util.logging.Level;

/**
 * Represents a point in three-dimensional Cartesian space.

 * This class encapsulates the coordinates (x, y, z) of a point and provides
 * utility methods for common geometric operations such as distance calculation,
 * rotation, and vector operations.

 * Design Patterns and Principles:

 * 1. VALUE OBJECT PATTERN: This class represents an immutable value object where
 *    the state (x, y, z coordinates) is set at construction and cannot be modified.
 *    Operations that would "change" the point instead return new Point3D instances.
 *    This ensures thread-safety and predictable behavior.

 * 2. ENCAPSULATION: The coordinates are stored as private final fields, accessible
 *    only through getter methods. This prevents external modification and maintains
 *    data integrity.

 * 3. FLUENT INTERFACE (partial): Methods return Point3D objects, allowing for
 *    method chaining in some contexts.

 * Algorithmic Foundations:
 * <p>
 * - Distance calculations use the Euclidean distance formula (Pythagorean theorem
 *   in 3D space), demonstrating fundamental geometric algorithms with O(1) complexity.
 * <p>
 * - Rotation operations use rotation matrices from linear algebra, showing how
 *   matrix transformations apply to spatial data structures.
 * <p>
 * - The class serves as a building block for more complex spatial data structures
 *   like octrees, k-d trees, or point clouds used in computational geometry.
 * 
 * @author Your Name
 * @version 1.0
 */
public class Point3D {
    
    private static final Logger logger = Logger.getLogger(Point3D.class.getName());


    private final double x;
    private final double y;
    private final double z;
    
    /**
     * Constructs a new Point3D with the specified coordinates.
     * <p>
     * This constructor initializes the point's position in 3D space. The immutability
     * of this object is guaranteed by making all fields final, which is a key principle
     * in creating thread-safe data structures.
     * 
     * @param x the x-coordinate
     * @param y the y-coordinate
     * @param z the z-coordinate
     */
    public Point3D(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
        logger.log(Level.INFO, "Created new Point3D at ({0}, {1}, {2})", 
                   new Object[]{x, y, z});
    }
    
    /**
     * Returns the x-coordinate of this point.
     * <p>
     * Getter methods are fundamental to encapsulation, providing controlled read-only
     * access to internal state. This ensures that the coordinate cannot be modified
     * after construction, maintaining the immutability contract.
     * 
     * @return the x-coordinate
     */
    public double getX() {
        return x;
    }
    
    /**
     * Returns the y-coordinate of this point.
     * 
     * @return the y-coordinate
     */
    public double getY() {
        return y;
    }
    
    /**
     * Returns the z-coordinate of this point.
     * 
     * @return the z-coordinate
     */
    public double getZ() {
        return z;
    }
    
    /**
     * Calculates the Euclidean distance from this point to another point.
     * <p>
     * This method implements the 3D Euclidean distance formula:
     * distance = sqrt((x2-x1)² + (y2-y1)² + (z2-z1)²)
     * <p>
     * The algorithm has O(1) time complexity and demonstrates fundamental geometric
     * calculations used in spatial indexing, nearest-neighbor searches, and clustering
     * algorithms.
     * <p>
     * Error handling: Returns Double.NaN if the other point is null, allowing the
     * caller to handle invalid input gracefully.
     * 
     * @param other the point to calculate distance to
     * @return the Euclidean distance between the two points, or Double.NaN if other is null
     */
    public double distanceTo(Point3D other) {
        if (other == null) {
            logger.log(Level.WARNING, "Attempted to calculate distance to null point");
            return Double.NaN;
        }
        
        double dx = this.x - other.x;
        double dy = this.y - other.y;
        double dz = this.z - other.z;
        double distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
        
        logger.log(Level.INFO, "Distance from ({0}, {1}, {2}) to ({3}, {4}, {5}) = {6}",
                   new Object[]{this.x, this.y, this.z, other.x, other.y, other.z, distance});
        
        return distance;
    }
    
    /**
     * Calculates the distance from this point to the origin (0, 0, 0).
     * <p>
     * This is a specialized case of distance calculation that computes the magnitude
     * or length of the position vector. It's commonly used in normalization operations
     * and when determining if points lie within a certain radius of the origin.
     * 
     * @return the distance from this point to the origin
     */
    public double distanceFromOrigin() {
        double distance = Math.sqrt(x * x + y * y + z * z);
        logger.log(Level.INFO, "Distance from origin for point ({0}, {1}, {2}) = {3}",
                   new Object[]{x, y, z, distance});
        return distance;
    }
    
    /**
     * Rotates this point around the X-axis by the specified angle.
     * <p>
     * This method applies a rotation matrix transformation:
     * [1    0         0    ] [x]
     * [0  cos(θ)  -sin(θ) ] [y]
     * [0  sin(θ)   cos(θ) ] [z]
     * <p>
     * The rotation follows the right-hand rule: positive angles rotate counterclockwise
     * when looking from positive X toward the origin.
     * <p>
     * Returns a new Point3D instance, maintaining immutability. This is a key principle
     * in functional programming and helps prevent bugs in concurrent environments.
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new Point3D representing the rotated point
     */
    public Point3D rotateX(double angleRadians) {
        double cos = Math.cos(angleRadians);
        double sin = Math.sin(angleRadians);
        
        double newY = y * cos - z * sin;
        double newZ = y * sin + z * cos;
        
        logger.log(Level.INFO, "Rotated point ({0}, {1}, {2}) around X-axis by {3} radians to ({4}, {5}, {6})",
                   new Object[]{x, y, z, angleRadians, x, newY, newZ});
        
        return new Point3D(x, newY, newZ);
    }
    
    /**
     * Rotates this point around the Y-axis by the specified angle.
     * <p>
     * This method applies a rotation matrix transformation:
     * [ cos(θ)  0  sin(θ)] [x]
     * [   0     1    0   ] [y]
     * [-sin(θ)  0  cos(θ)] [z]
     * <p>
     * The rotation follows the right-hand rule: positive angles rotate counterclockwise
     * when looking from positive Y toward the origin.
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new Point3D representing the rotated point
     */
    public Point3D rotateY(double angleRadians) {
        double cos = Math.cos(angleRadians);
        double sin = Math.sin(angleRadians);
        
        double newX = x * cos + z * sin;
        double newZ = -x * sin + z * cos;
        
        logger.log(Level.INFO, "Rotated point ({0}, {1}, {2}) around Y-axis by {3} radians to ({4}, {5}, {6})",
                   new Object[]{x, y, z, angleRadians, newX, y, newZ});
        
        return new Point3D(newX, y, newZ);
    }
    
    /**
     * Rotates this point around the Z-axis by the specified angle.
     * <p>
     * This method applies a rotation matrix transformation:
     * [cos(θ)  -sin(θ)  0] [x]
     * [sin(θ)   cos(θ)  0] [y]
     * [  0        0     1] [z]
     * <p>
     * The rotation follows the right-hand rule: positive angles rotate counterclockwise
     * when looking from positive Z toward the origin.
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new Point3D representing the rotated point
     */
    public Point3D rotateZ(double angleRadians) {
        double cos = Math.cos(angleRadians);
        double sin = Math.sin(angleRadians);
        
        double newX = x * cos - y * sin;
        double newY = x * sin + y * cos;
        
        logger.log(Level.INFO, "Rotated point ({0}, {1}, {2}) around Z-axis by {3} radians to ({4}, {5}, {6})",
                   new Object[]{x, y, z, angleRadians, newX, newY, z});
        
        return new Point3D(newX, newY, z);
    }
    
    /**
     * Calculates the midpoint between this point and another point.
     * <p>
     * The midpoint is calculated as the average of corresponding coordinates:
     * midpoint = ((x1+x2)/2, (y1+y2)/2, (z1+z2)/2)
     * <p>
     * This is commonly used in divide-and-conquer algorithms, binary space partitioning,
     * and interpolation techniques.
     * 
     * @param other the other point
     * @return a new Point3D representing the midpoint, or null if other is null
     */
    public Point3D midpoint(Point3D other) {
        if (other == null) {
            logger.log(Level.WARNING, "Attempted to calculate midpoint with null point");
            return null;
        }
        
        double midX = (this.x + other.x) / 2.0;
        double midY = (this.y + other.y) / 2.0;
        double midZ = (this.z + other.z) / 2.0;
        
        logger.log(Level.INFO, "Calculated midpoint between ({0}, {1}, {2}) and ({3}, {4}, {5})",
                   new Object[]{this.x, this.y, this.z, other.x, other.y, other.z});
        
        return new Point3D(midX, midY, midZ);
    }
    
    /**
     * Translates this point by the specified offset values.
     * <p>
     * Translation is a fundamental affine transformation that moves a point by adding
     * offset values to each coordinate. This operation has O(1) complexity and is
     * frequently used in animation, physics simulations, and coordinate system conversions.
     * 
     * @param dx the offset in the x direction
     * @param dy the offset in the y direction
     * @param dz the offset in the z direction
     * @return a new Point3D representing the translated point
     */
    public Point3D translate(double dx, double dy, double dz) {
        logger.log(Level.INFO, "Translating point ({0}, {1}, {2}) by ({3}, {4}, {5})",
                   new Object[]{x, y, z, dx, dy, dz});
        return new Point3D(x + dx, y + dy, z + dz);
    }
    
    /**
     * Scales this point by the specified factor relative to the origin.
     * <p>
     * Scaling multiplies each coordinate by the scale factor. This is useful for
     * zooming operations, normalization, and coordinate system transformations.
     * <p>
     * Warning: A scale factor of 0 will collapse the point to the origin.
     * 
     * @param scale the scaling factor
     * @return a new Point3D representing the scaled point
     */
    public Point3D scale(double scale) {
        if (scale == 0) {
            logger.log(Level.WARNING, "Scaling point ({0}, {1}, {2}) by 0 - point will collapse to origin",
                       new Object[]{x, y, z});
        } else {
            logger.log(Level.INFO, "Scaling point ({0}, {1}, {2}) by factor {3}",
                       new Object[]{x, y, z, scale});
        }
        return new Point3D(x * scale, y * scale, z * scale);
    }
    
    /**
     * Calculates the dot product of this point (treated as a vector) with another point.
     * <p>
     * The dot product is defined as: x1*x2 + y1*y2 + z1*z2
     * <p>
     * This operation is fundamental in computational geometry and has applications in:
     * - Determining the angle between vectors
     * - Projecting one vector onto another
     * - Testing for orthogonality (dot product = 0)
     * 
     * @param other the other point/vector
     * @return the dot product, or Double.NaN if other is null
     */
    public double dotProduct(Point3D other) {
        if (other == null) {
            logger.log(Level.WARNING, "Attempted to calculate dot product with null point");
            return Double.NaN;
        }
        
        double result = this.x * other.x + this.y * other.y + this.z * other.z;
        logger.log(Level.INFO, "Dot product of ({0}, {1}, {2}) and ({3}, {4}, {5}) = {6}",
                   new Object[]{this.x, this.y, this.z, other.x, other.y, other.z, result});
        
        return result;
    }
    
    /**
     * Calculates the cross product of this point (treated as a vector) with another point.
     * <p>
     * The cross product produces a vector perpendicular to both input vectors:
     * result = (y1*z2 - z1*y2, z1*x2 - x1*z2, x1*y2 - y1*x2)
     * <p>
     * This operation is essential for:
     * - Calculating surface normals in 3D graphics
     * - Determining the area of parallelograms
     * - Finding perpendicular vectors
     * 
     * @param other the other point/vector
     * @return a new Point3D representing the cross product, or null if other is null
     */
    public Point3D crossProduct(Point3D other) {
        if (other == null) {
            logger.log(Level.WARNING, "Attempted to calculate cross product with null point");
            return null;
        }
        
        double newX = this.y * other.z - this.z * other.y;
        double newY = this.z * other.x - this.x * other.z;
        double newZ = this.x * other.y - this.y * other.x;
        
        logger.log(Level.INFO, "Cross product of ({0}, {1}, {2}) and ({3}, {4}, {5})",
                   new Object[]{this.x, this.y, this.z, other.x, other.y, other.z});
        
        return new Point3D(newX, newY, newZ);
    }
    
    /**
     * Normalizes this point (treated as a vector) to unit length.
     * <p>
     * Normalization scales the vector so that its magnitude equals 1 while preserving
     * its direction. This is crucial in graphics, physics, and many geometric algorithms.
     * <p>
     * The normalized vector is calculated by dividing each component by the vector's magnitude.
     * <p>
     * Error handling: If the point is at the origin (magnitude = 0), normalization is
     * mathematically undefined. This method returns null and logs a severe error.
     * 
     * @return a new Point3D representing the normalized vector, or null if magnitude is 0
     */
    public Point3D normalize() {
        double magnitude = distanceFromOrigin();
        
        if (magnitude == 0) {
            logger.log(Level.SEVERE, "Cannot normalize zero vector at origin");
            return null;
        }
        
        logger.log(Level.INFO, "Normalizing vector ({0}, {1}, {2}) with magnitude {3}",
                   new Object[]{x, y, z, magnitude});
        
        return new Point3D(x / magnitude, y / magnitude, z / magnitude);
    }
    
    /**
     * Checks if this point is equal to another object.
     * <p>
     * Two points are considered equal if their coordinates are equal within a small
     * epsilon tolerance (1e-10) to account for floating-point precision issues.
     * <p>
     * This demonstrates proper implementation of the equals() method contract, which
     * is essential for using objects in hash-based collections and for correctness
     * in data structure operations.
     * 
     * @param obj the object to compare with
     * @return true if the points are equal within tolerance, false otherwise
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        
        Point3D other = (Point3D) obj;
        double epsilon = 1e-10;
        
        boolean equal = Math.abs(this.x - other.x) < epsilon &&
                       Math.abs(this.y - other.y) < epsilon &&
                       Math.abs(this.z - other.z) < epsilon;
        
        if (!equal) {
            logger.log(Level.INFO, "Point equality check: ({0}, {1}, {2}) != ({3}, {4}, {5})",
                       new Object[]{this.x, this.y, this.z, other.x, other.y, other.z});
        }
        
        return equal;
    }
    
    /**
     * Returns a hash code for this point.
     * <p>
     * This method must be overridden whenever equals() is overridden to maintain
     * the general contract for hashCode, which states that equal objects must have
     * equal hash codes. This is critical for correct behavior in hash-based collections
     * like HashMap and HashSet.
     * <p>
     * Uses the Objects.hash() utility method for a robust hash code calculation.
     * 
     * @return a hash code value for this point
     */
    @Override
    public int hashCode() {
        return java.util.Objects.hash(x, y, z);
    }
    
    /**
     * Returns a string representation of this point.
     * <p>
     * The string format "Point3D(x, y, z)" provides clear, readable output for
     * debugging and logging purposes. This follows the common convention for
     * toString() implementations.
     * 
     * @return a string representation of this point
     */
    @Override
    public String toString() {
        return String.format("Point3D(%.2f, %.2f, %.2f)", x, y, z);
    }

    public Point3D subtract(Point3D other) {
        if (other == null) {
            logger.log(Level.WARNING, "Attempted to subtract null point");
            return null; // or throw IllegalArgumentException
        }

        return new Point3D(this.x - other.x, this.y - other.y, this.z - other.z);
    }

        public double distanceSquaredTo (Point3D other){
            if (other == null) {
                logger.log(Level.WARNING, "Attempted to calculate squared distance to null point");
                return Double.NaN;
            }

            double dx = this.x - other.x;
            double dy = this.y - other.y;
            double dz = this.z - other.z;

            return dx * dx + dy * dy + dz * dz;
        }
    }
