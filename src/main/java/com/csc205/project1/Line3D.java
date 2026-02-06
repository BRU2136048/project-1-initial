package com.csc205.project1;

import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

    /**
     * Represents a line segment in 3D space, defined by two endpoints.
     *
     * <p>Design note:
     * - This class models a <b>line segment</b> (finite length) because it has a meaningful {@link #length()}.
     * - It can also provide computations for the corresponding <b>infinite line</b> that extends through the endpoints.
     *
     * <p>Dependencies:
     * - This class expects a previously-created {@code Point3D} class (from your earlier work).
     *   Typical assumed API:
     *     - double getX(), getY(), getZ()
     *     - double distanceTo(Point3D other)
     */
    public final class Line3D {

        private static final Logger logger = Logger.getLogger(Line3D.class.getName());

        private final Point3D start;
        private final Point3D end;

        /**
         * Creates a new 3D line segment from {@code start} to {@code end}.
         *
         * <p>This constructor validates that both endpoints are provided. It does not forbid
         * a "degenerate" segment (start == end), but methods will log warnings where it matters.
         *
         * @param start the starting endpoint of the segment
         * @param end   the ending endpoint of the segment
         */
        public Line3D(Point3D start, Point3D end) {
            this.start = Objects.requireNonNull(start, "start must not be null");
            this.end = Objects.requireNonNull(end, "end must not be null");

            logger.log(Level.INFO, () -> "Created Line3D with start=" + start + ", end=" + end);

            if (start.equals(end)) {
                logger.log(Level.WARNING, "Line3D is degenerate (start == end). Length is 0 and direction is undefined.");
            }
        }

        /**
         * Returns the start endpoint of this line segment.
         *
         * <p>Exposing the endpoint via a getter supports clean encapsulation: callers can
         * read the state needed for algorithms (e.g., intersection checks) without mutating
         * this object.
         *
         * @return the start endpoint
         */
        public Point3D getStart() {
            logger.log(Level.INFO, "getStart() called.");
            return start;
        }

        /**
         * Returns the end endpoint of this line segment.
         *
         * <p>Like {@link #getStart()}, this accessor enables algorithmic operations on the segment
         * while keeping the segment immutable.
         *
         * @return the end endpoint
         */
        public Point3D getEnd() {
            logger.log(Level.INFO, "getEnd() called.");
            return end;
        }

        /**
         * Computes the length of this line segment.
         *
         * <p>This method uses the distance function on {@code Point3D} (typically Euclidean distance).
         * The runtime is constant O(1) because it performs a fixed number of arithmetic operations.
         *
         * @return the segment length (0.0 for a degenerate segment)
         */
        public double length() {
            logger.log(Level.INFO, "Computing segment length.");
            try {
                double len = start.distanceTo(end);
                logger.log(Level.INFO, () -> "Segment length computed: " + len);
                return len;
            } catch (RuntimeException ex) {
                logger.log(Level.SEVERE, "Failed computing length due to unexpected error.", ex);
                throw ex;
            }
        }

        /**
         * Returns the midpoint of this line segment.
         *
         * <p>The midpoint is often useful in geometry algorithms and as a stable representative
         * of the segment for heuristics (e.g., spatial partitioning).
         *
         * @return the midpoint between start and end
         */
        public Point3D midpoint() {
            logger.log(Level.INFO, "Computing midpoint.");
            try {
                double mx = (start.getX() + end.getX()) / 2.0;
                double my = (start.getY() + end.getY()) / 2.0;
                double mz = (start.getZ() + end.getZ()) / 2.0;

                Point3D mid = new Point3D(mx, my, mz);
                logger.log(Level.INFO, () -> "Midpoint computed: " + mid);
                return mid;
            } catch (RuntimeException ex) {
                logger.log(Level.SEVERE, "Failed computing midpoint due to unexpected error.", ex);
                throw ex;
            }
        }

        /**
         * Computes the direction vector of the infinite line passing through this segment.
         *
         * <p>This vector points from {@code start} to {@code end}. Many 3D algorithms (distance between lines,
         * projections, intersections) rely on direction vectors plus dot/cross products.
         *
         * <p>If the segment is degenerate, the direction is the zero vector and a warning is logged.
         *
         * @return a {@link Vector3} representing (end - start)
         */
        public Vector3 direction() {
            logger.log(Level.INFO, "Computing direction vector (end - start).");
            Vector3 d = Vector3.fromPoints(start, end);
            if (d.isZero()) {
                logger.log(Level.WARNING, "Direction is zero vector (degenerate segment).");
            } else {
                logger.log(Level.INFO, () -> "Direction vector computed: " + d);
            }
            return d;
        }

        /**
         * Checks whether this line's infinite extension is parallel to another line's infinite extension.
         *
         * <p>Two lines are parallel when their direction vectors are scalar multiples. In 3D,
         * we can detect this by checking whether the cross product of the direction vectors is near zero.
         *
         * @param other the other line
         * @return true if the infinite lines are parallel (within a small tolerance), false otherwise
         */
        public boolean isParallelTo(Line3D other) {
            Objects.requireNonNull(other, "other must not be null");
            logger.log(Level.INFO, "Checking parallelism with another line.");

            Vector3 d1 = this.direction();
            Vector3 d2 = other.direction();

            if (d1.isZero() || d2.isZero()) {
                logger.log(Level.WARNING, "Parallel check involves degenerate line(s). Treating zero-direction as parallel only if both are zero.");
                return d1.isZero() && d2.isZero();
            }

            boolean parallel = d1.cross(d2).isNearZero(Vector3.EPS);
            logger.log(Level.INFO, () -> "isParallelTo result: " + parallel);
            return parallel;
        }

        /**
         * Computes the shortest distance between the <b>infinite lines</b> defined by this segment and {@code other}.
         *
         * <p>This is the classic "distance between skew lines" computation:
         * <ul>
         *   <li>If lines are not parallel, distance = |(p2 - p1) · (d1 × d2)| / |d1 × d2|</li>
         *   <li>If lines are parallel, distance = distance from p2 to line1</li>
         * </ul>
         *
         * <p>Algorithmic note:
         * This is O(1) time and O(1) space: it uses a fixed set of vector operations.
         *
         * @param other another 3D line segment whose infinite extension is considered
         * @return the shortest distance between the two infinite lines
         */
        public double shortestDistanceBetweenInfiniteLines(Line3D other) {
            Objects.requireNonNull(other, "other must not be null");
            logger.log(Level.INFO, "Computing shortest distance between infinite lines.");

            Vector3 d1 = this.direction();
            Vector3 d2 = other.direction();

            // Handle degenerate cases early.
            if (d1.isZero() && d2.isZero()) {
                logger.log(Level.WARNING, "Both lines are degenerate points; returning point-to-point distance.");
                return this.start.distanceTo(other.start);
            }
            if (d1.isZero()) {
                logger.log(Level.WARNING, "This line is degenerate; returning distance from this point to other infinite line.");
                return distanceFromPointToInfiniteLine(this.start, other);
            }
            if (d2.isZero()) {
                logger.log(Level.WARNING, "Other line is degenerate; returning distance from other point to this infinite line.");
                return distanceFromPointToInfiniteLine(other.start, this);
            }

            Vector3 n = d1.cross(d2);                 // Normal to both directions.
            double nMag = n.magnitude();

            if (nMag < Vector3.EPS) {
                // Parallel lines: distance from other.start to this infinite line.
                logger.log(Level.INFO, "Lines are parallel (cross magnitude near zero). Using point-to-line distance.");
                double dist = distanceFromPointToInfiniteLine(other.start, this);
                logger.log(Level.INFO, () -> "Parallel-line distance computed: " + dist);
                return dist;
            }

            Vector3 p21 = Vector3.fromPoints(this.start, other.start); // (other.start - this.start)
            double dist = Math.abs(p21.dot(n)) / nMag;

            logger.log(Level.INFO, () -> "Shortest distance between infinite lines computed: " + dist);
            return dist;
        }

        /**
         * Computes the shortest distance between the <b>line segments</b> represented by this line and {@code other}.
         *
         * <p>This method solves the closest points between two segments in 3D using a standard
         * clamp-to-[0,1] approach. It is robust for parallel and non-parallel configurations.
         *
         * <p>Practical note:
         * Segment-to-segment distance is commonly used in collision detection, CAD, and simulation.
         *
         * @param other another line segment
         * @return the shortest distance between the two segments
         */
        public double shortestDistanceBetweenSegments(Line3D other) {
            Objects.requireNonNull(other, "other must not be null");
            logger.log(Level.INFO, "Computing shortest distance between line segments.");

            Vector3   u = Vector3.fromPoints(this.start, this.end);
            Vector3   v = Vector3.fromPoints(other.start, other.end);
            Vector3   w = Vector3.fromPoints(other.start, this.start); // (this.start - other.start)

            double a = u.dot(u); // |u|^2
            double b = u.dot(v);
            double c = v.dot(v); // |v|^2
            double d = u.dot(w);
            double e = v.dot(w);

            // Degenerate segments: treat as point-segment or point-point.
            if (a < Vector3.EPS && c < Vector3.EPS) {
                logger.log(Level.WARNING, "Both segments degenerate into points; returning point-to-point distance.");
                return this.start.distanceTo(other.start);
            }
            if (a < Vector3.EPS) {
                logger.log(Level.WARNING, "This segment degenerates into a point; returning point-to-segment distance.");
                return distanceFromPointToSegment(this.start, other);
            }
            if (c < Vector3.EPS) {
                logger.log(Level.WARNING, "Other segment degenerates into a point; returning point-to-segment distance.");
                return distanceFromPointToSegment(other.start, this);
            }

            double denom = a * c - b * b; // Always >= 0

            double sN, sD = denom;
            double tN, tD = denom;

            // If denom is 0, lines are (almost) parallel.
            if (denom < Vector3.EPS) {
                logger.log(Level.INFO, "Segments are nearly parallel. Falling back to parallel handling.");
                sN = 0.0;
                sD = 1.0;
                tN = e;
                tD = c;
            } else {
                sN = (b * e - c * d);
                tN = (a * e - b * d);

                // Clamp s to [0,1]
                if (sN < 0.0) {
                    sN = 0.0;
                    tN = e;
                    tD = c;
                } else if (sN > sD) {
                    sN = sD;
                    tN = e + b;
                    tD = c;
                }
            }

            // Clamp t to [0,1]
            if (tN < 0.0) {
                tN = 0.0;
                // Recompute s for this edge case
                if (-d < 0.0) {
                    sN = 0.0;
                } else if (-d > a) {
                    sN = sD;
                } else {
                    sN = -d;
                    sD = a;
                }
            } else if (tN > tD) {
                tN = tD;
                // Recompute s for this edge case
                if ((-d + b) < 0.0) {
                    sN = 0.0;
                } else if ((-d + b) > a) {
                    sN = sD;
                } else {
                    sN = (-d + b);
                    sD = a;
                }
            }

            double s = (Math.abs(sN) < Vector3.EPS) ? 0.0 : (sN / sD);
            double t = (Math.abs(tN) < Vector3.EPS) ? 0.0 : (tN / tD);

            Vector3 closestDelta = w.add(u.scale(s)).subtract(v.scale(t));
            double dist = closestDelta.magnitude();

            logger.log(Level.INFO, () -> "Shortest distance between segments computed: " + dist);
            return dist;
        }

        /**
         * Computes the distance from {@code p} to the infinite line defined by {@code line}.
         *
         * <p>This is used as a building block for handling parallel lines and degenerate cases.
         *
         * @param p    a point
         * @param line a line whose infinite extension is considered
         * @return the shortest distance from the point to the infinite line
         */
        public static double distanceFromPointToInfiniteLine(Point3D p, Line3D line) {
            Objects.requireNonNull(p, "p must not be null");
            Objects.requireNonNull(line, "line must not be null");

            logger.log(Level.INFO, "Computing distance from point to infinite line.");

            Vector3 d = line.direction();
            if (d.isZero()) {
                logger.log(Level.WARNING, "Target line is degenerate; returning point-to-point distance.");
                return p.distanceTo(line.start);
            }

            Vector3 ap = Vector3.fromPoints(line.start, p); // (p - line.start)
            // Distance is |ap x d| / |d|
            double dist = ap.cross(d).magnitude() / d.magnitude();

            logger.log(Level.INFO, () -> "Point-to-infinite-line distance computed: " + dist);
            return dist;
        }

        /**
         * Computes the distance from {@code p} to the finite segment {@code seg}.
         *
         * <p>This is a standard projection-and-clamp algorithm:
         * <ol>
         *   <li>Project (p - start) onto the segment direction.</li>
         *   <li>Clamp the parameter to [0,1] to stay on the segment.</li>
         *   <li>Measure distance to the resulting closest point.</li>
         * </ol>
         *
         * @param p   a point
         * @param seg a segment
         * @return the shortest distance from the point to the segment
         */
        public static double distanceFromPointToSegment(Point3D p, Line3D seg) {
            Objects.requireNonNull(p, "p must not be null");
            Objects.requireNonNull(seg, "seg must not be null");

            logger.log(Level.INFO, "Computing distance from point to segment.");

            Vector3 ab = Vector3.fromPoints(seg.start, seg.end);
            double ab2 = ab.dot(ab);

            if (ab2 < Vector3.EPS) {
                logger.log(Level.WARNING, "Segment is degenerate; returning point-to-point distance.");
                return p.distanceTo(seg.start);
            }

            Vector3 ap = Vector3.fromPoints(seg.start, p);
            double t = ap.dot(ab) / ab2;
            t = clamp(t, 0.0, 1.0);

            Point3D closest = new Point3D(
                    seg.start.getX() + ab.x * t,
                    seg.start.getY() + ab.y * t,
                    seg.start.getZ() + ab.z * t
            );

            double dist = p.distanceTo(closest);
            logger.log(Level.INFO, () -> "Point-to-segment distance computed: " + dist + " (closest=" + closest + ")");
            return dist;
        }

        /**
         * Provides a concise textual representation useful for logging and debugging.
         *
         * @return a string describing the segment endpoints
         */
        @Override
        public String toString() {
            return "Line3D{start=" + start + ", end=" + end + '}';
        }

        /**
         * Two Line3D instances are equal if their endpoints are equal in the same order.
         *
         * <p>Note:
         * If you want undirected equality (A->B equals B->A), you can adjust this implementation.
         *
         * @param o the other object
         * @return true if equal, false otherwise
         */
        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof Line3D)) return false;
            Line3D line3D = (Line3D) o;
            return start.equals(line3D.start) && end.equals(line3D.end);
        }

        /**
         * Hash code consistent with {@link #equals(Object)}.
         *
         * @return hash code for this instance
         */
        @Override
        public int hashCode() {
            return Objects.hash(start, end);
        }

        // ----------------------------
        // Helper types and utilities
        // ----------------------------

        /**
         * Minimal 3D vector helper used internally by Line3D.
         *
         * <p>This intentionally keeps vector math close to the algorithms that use it, making the
         * data flow easy to follow (a common goal in Spring getting-started style examples).
         */
        public static final class Vector3 {
            // A small epsilon for floating-point comparisons in geometry.
            public static final double EPS = 1e-10;

            public final double x;
            public final double y;
            public final double z;

            public Vector3(double x, double y, double z) {
                this.x = x;
                this.y = y;
                this.z = z;
            }

            /**
             * Builds a vector from {@code a} to {@code b}, i.e., (b - a).
             *
             * @param a start point
             * @param b end point
             * @return vector (b - a)
             */
            public static Vector3 fromPoints(Point3D a, Point3D b) {
                return new Vector3(b.getX() - a.getX(), b.getY() - a.getY(), b.getZ() - a.getZ());
            }

            /**
             * Dot product between this vector and {@code other}.
             *
             * <p>The dot product is a core operation for projections and angle relationships.
             *
             * @param other another vector
             * @return scalar dot product
             */
            public double dot(Vector3 other) {
                return this.x * other.x + this.y * other.y + this.z * other.z;
            }

            /**
             * Cross product between this vector and {@code other}.
             *
             * <p>The cross product produces a vector orthogonal to both inputs. It is essential
             * for computing distances between skew lines and detecting parallel direction vectors.
             *
             * @param other another vector
             * @return cross product vector
             */
            public Vector3 cross(Vector3 other) {
                return new Vector3(
                        this.y * other.z - this.z * other.y,
                        this.z * other.x - this.x * other.z,
                        this.x * other.y - this.y * other.x
                );
            }

            /**
             * Vector magnitude (length).
             *
             * @return Euclidean norm
             */
            public double magnitude() {
                return Math.sqrt(this.dot(this));
            }

            /**
             * Scales this vector by {@code s}.
             *
             * @param s scalar value
             * @return scaled vector
             */
            public Vector3 scale(double s) {
                return new Vector3(this.x * s, this.y * s, this.z * s);
            }

            /**
             * Adds {@code other} to this vector.
             *
             * @param other another vector
             * @return sum vector
             */
            public Vector3 add(Vector3 other) {
                return new Vector3(this.x + other.x, this.y + other.y, this.z + other.z);
            }

            /**
             * Subtracts {@code other} from this vector.
             *
             * @param other another vector
             * @return difference vector
             */
            public Vector3 subtract(Vector3 other) {
                return new Vector3(this.x - other.x, this.y - other.y, this.z - other.z);
            }

            /**
             * Checks whether this is (approximately) the zero vector.
             *
             * @return true if vector magnitude is near zero
             */
            public boolean isZero() {
                return isNearZero(EPS);
            }

            /**
             * Checks whether this vector is near zero using a supplied tolerance.
             *
             * @param eps tolerance value
             * @return true if magnitude < eps
             */
            public boolean isNearZero(double eps) {
                return Math.abs(x) < eps && Math.abs(y) < eps && Math.abs(z) < eps;
            }

            @Override
            public String toString() {
                return "Vector3{" + "x=" + x + ", y=" + y + ", z=" + z + '}';
            }
        }

        private static double clamp(double value, double min, double max) {
            if (value < min) return min;
            if (value > max) return max;
            return value;
        }

        // ---------------------------------------------
        // Object-oriented design patterns & DSA notes
        // ---------------------------------------------

        /**
         * Documents the design patterns and how they connect to foundational DSA principles.
         *
         * <p>You can print or expose this string in your coursework, README, or Javadoc.
         *
         * @return explanation of design choices and their DSA relevance
         */
        public static String designPatternsAndDsaRationale() {
            return ""
                    + "Design patterns used:\n"
                    + "1) Immutable Value Object (also aligns with 'Value Object' concept):\n"
                    + "   - Line3D is final, fields are private final, and there are no setters.\n"
                    + "   - Benefits: thread-safety by design, easier reasoning, fewer state bugs.\n"
                    + "   - DSA connection: many geometric algorithms assume inputs do not mutate\n"
                    + "     mid-computation. Immutability keeps algorithm invariants stable.\n\n"
                    + "2) Composition over Inheritance:\n"
                    + "   - Line3D is composed of two Point3D instances.\n"
                    + "   - Benefits: clearer data modeling and reuse of Point3D operations.\n"
                    + "   - DSA connection: composition makes it straightforward to build higher-level\n"
                    + "     data structures (e.g., polylines, meshes, graphs of edges) from points/edges.\n\n"
                    + "3) Helper/Utility Type (Vector3) scoped to the domain:\n"
                    + "   - Vector3 encapsulates dot/cross/projection primitives.\n"
                    + "   - Benefits: separates low-level operations from higher-level intent.\n"
                    + "   - DSA connection: dot/cross products are the primitive operations that power\n"
                    + "     higher-level algorithms (closest approach, collision checks, spatial indexing).\n\n"
                    + "4) Defensive Programming (Guard Clauses) + Logging as Observability:\n"
                    + "   - Null checks and degenerate-case handling are explicit.\n"
                    + "   - INFO/WARNING/SEVERE logs communicate algorithm flow and failure points.\n"
                    + "   - DSA connection: edge-case correctness is a major part of algorithmic rigor.\n\n"
                    + "Foundational principles demonstrated:\n"
                    + "- Encapsulation: internal math details are hidden behind expressive methods.\n"
                    + "- Abstraction: callers think in terms of 'distance between lines' not cross products.\n"
                    + "- Algorithmic efficiency: all core computations are O(1) time and O(1) space.\n"
                    + "- Numerical robustness: EPS tolerances acknowledge floating-point behavior.\n";
        }
    }


