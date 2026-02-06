package com.csc205.project1;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * A 3D cube suitable for graphics / simulation code.
 *
 * <p>Design goals:
 * <ul>
 *   <li>Represent the cube via its 8 vertices (a compact, graphics-friendly structure).</li>
 *   <li>Expose common transforms (translate / rotate) as methods that mutate this cube.</li>
 *   <li>Provide geometry queries (volume, surface area, total edge length, edges, vertices).</li>
 * </ul>
 *
 * <p>Dependencies:
 * <ul>
 *   <li>{@code Point3D}: a 3D point value type (x,y,z)</li>
 *   <li>{@code Line3D}: a line segment value type (start,end) with a length() method</li>
 * </ul>
 */
public class Cube3D {

    private static final Logger log = Logger.getLogger(Cube3D.class.getName());

    /**
     * Axis selector for convenience rotations.
     * This avoids "magic strings" and keeps rotation calls type-safe.
     */
    public enum Axis {
        X, Y, Z
    }

    // --- Core state ---
    private Point3D center;
    private double sideLength;

    // Stored vertices (in a stable order) to support fast transforms and edge construction.
    // Vertex order is defined by (xSign, ySign, zSign) for +/- half side length.
    private final Point3D[] vertices = new Point3D[8];

    // Optional cached edges; rebuilt on demand if null.
    private List<Line3D> cachedEdges;

    /**
     * Create a cube from a center point and side length.
     *
     * <p>In a graphics program, "center + size" is a common representation because it’s easy
     * to transform and easy to compute bounding volumes for. This constructor expands that
     * compact representation into 8 vertices.</p>
     *
     * @param center the cube center
     * @param sideLength the edge length; must be positive
     */
    public Cube3D(Point3D center, double sideLength) {
        this.center = Objects.requireNonNull(center, "center must not be null");
        if (sideLength <= 0.0) {
            log.severe("Attempted to create Cube3D with non-positive sideLength: " + sideLength);
            throw new IllegalArgumentException("sideLength must be > 0");
        }
        this.sideLength = sideLength;

        log.info(() -> "Creating Cube3D center=" + center + ", sideLength=" + sideLength);
        rebuildVerticesFromCenter();
    }

    /**
     * Factory method that builds a cube from its minimum and maximum corners.
     *
     * <p>This is useful when your pipeline produces axis-aligned bounds (AABB) and you need
     * to convert them into a cube instance. The method validates that the provided box is
     * actually a cube (equal extents in X/Y/Z) and then creates the Cube3D.</p>
     *
     * @param minCorner the minimum (x,y,z)
     * @param maxCorner the maximum (x,y,z)
     * @return a new Cube3D representing that cube-shaped volume
     */
    public static Cube3D fromMinMax(Point3D minCorner, Point3D maxCorner) {
        Objects.requireNonNull(minCorner, "minCorner must not be null");
        Objects.requireNonNull(maxCorner, "maxCorner must not be null");

        double dx = maxCorner.getX() - minCorner.getX();
        double dy = maxCorner.getY() - minCorner.getY();
        double dz = maxCorner.getZ() - minCorner.getZ();

        // Floating point tolerant equality check
        double eps = 1e-9;
        if (dx <= 0 || dy <= 0 || dz <= 0 || (Math.abs(dx - dy) > eps) || (Math.abs(dy - dz) > eps)) {
            Logger.getLogger(Cube3D.class.getName()).severe(
                    "fromMinMax failed: not a valid cube. dx=" + dx + ", dy=" + dy + ", dz=" + dz);
            throw new IllegalArgumentException("min/max corners must describe a cube with positive equal extents");
        }

        Point3D center = new Point3D(
                (minCorner.getX() + maxCorner.getX()) / 2.0,
                (minCorner.getY() + maxCorner.getY()) / 2.0,
                (minCorner.getZ() + maxCorner.getZ()) / 2.0
        );

        Logger.getLogger(Cube3D.class.getName()).info(
                () -> "fromMinMax building Cube3D center=" + center + ", sideLength=" + dx);

        return new Cube3D(center, dx);
    }

    /**
     * Return the cube center.
     *
     * <p>Graphics engines often use a center point as the pivot for rotations and as the anchor
     * for spatial data structures (octrees, BVHs). This getter is intentionally cheap (O(1)).</p>
     */
    public Point3D getCenter() {
        return center;
    }

    /**
     * Return the cube edge length.
     *
     * <p>Even after rotations and translations, the cube remains rigid, so the side length is
     * stable unless explicitly changed.</p>
     */
    public double getSideLength() {
        return sideLength;
    }

    /**
     * Return an immutable snapshot list of the cube vertices (8 points).
     *
     * <p>Vertex arrays are the standard input for rendering pipelines. This method returns a
     * defensive copy so callers can’t accidentally corrupt the cube’s internal state.</p>
     */
    public List<Point3D> getVertices() {
        List<Point3D> out = new ArrayList<>(8);
        Collections.addAll(out, vertices);
        return Collections.unmodifiableList(out);
    }

    /**
     * Return the 12 cube edges as {@code Line3D} segments.
     *
     * <p>Edges are constructed from the vertex list using a fixed index map. The result is cached
     * for efficiency and invalidated whenever the cube is transformed.</p>
     */
    public List<Line3D> getEdges() {
        if (cachedEdges == null) {
            log.info("Rebuilding cached edges for Cube3D.");
            cachedEdges = Collections.unmodifiableList(buildEdges());
        }
        return cachedEdges;
    }

    /**
     * Translate the cube by a delta vector.
     *
     * <p>Translation is one of the most common transforms in a 3D program. This method shifts the
     * center and all vertices by (dx, dy, dz). The cube stays congruent; only its position changes.</p>
     *
     * @param dx delta X
     * @param dy delta Y
     * @param dz delta Z
     */
    public void translate(double dx, double dy, double dz) {
        log.info(() -> "Translating Cube3D by dx=" + dx + ", dy=" + dy + ", dz=" + dz);

        if (!isFinite(dx) || !isFinite(dy) || !isFinite(dz)) {
            log.severe("translate failed due to non-finite delta(s).");
            throw new IllegalArgumentException("Translation deltas must be finite numbers.");
        }

        center = new Point3D(center.getX() + dx, center.getY() + dy, center.getZ() + dz);

        for (int i = 0; i < vertices.length; i++) {
            Point3D v = vertices[i];
            vertices[i] = new Point3D(v.getX() + dx, v.getY() + dy, v.getZ() + dz);
        }

        invalidateCache();
    }

    /**
     * Rotate the cube around its center on a given axis.
     *
     * <p>This is the "game engine default" rotation: rotate around the cube’s own pivot (its center).
     * Under the hood, this is a matrix multiply applied to each vertex relative to the pivot.</p>
     *
     * @param axis the axis to rotate around
     * @param radians the rotation angle in radians
     */
    public void rotate(Axis axis, double radians) {
        Objects.requireNonNull(axis, "axis must not be null");
        rotateAroundPivot(axis, radians, this.center);
    }

    /**
     * Rotate the cube around an arbitrary pivot point on a given axis.
     *
     * <p>In 3D graphics, rotating around arbitrary pivots is essential (think: rotating a door around a hinge,
     * or orbiting an object around a point). This method supports that directly by subtracting the pivot,
     * applying rotation, then adding the pivot back.</p>
     *
     * @param axis the axis to rotate around
     * @param radians the rotation angle in radians
     * @param pivot the pivot point in world coordinates
     */
    public void rotateAroundPivot(Axis axis, double radians, Point3D pivot) {
        Objects.requireNonNull(axis, "axis must not be null");
        Objects.requireNonNull(pivot, "pivot must not be null");

        log.info(() -> "Rotating Cube3D around axis=" + axis + " radians=" + radians + " pivot=" + pivot);

        if (!isFinite(radians)) {
            log.severe("rotateAroundPivot failed due to non-finite angle.");
            throw new IllegalArgumentException("Rotation angle must be finite.");
        }

        double[][] R = rotationMatrix(axis, radians);

        // Rotate center as well (important when pivot != center).
        center = rotatePoint(center, pivot, R);

        for (int i = 0; i < vertices.length; i++) {
            vertices[i] = rotatePoint(vertices[i], pivot, R);
        }

        invalidateCache();
    }

    /**
     * Compute the cube volume.
     *
     * <p>Volume is a constant for rigid cubes and can be used for physics (mass from density),
     * collision heuristics, and sanity checks in geometry pipelines.</p>
     */
    public double volume() {
        double v = sideLength * sideLength * sideLength;
        log.info(() -> "Cube3D volume computed: " + v);
        return v;
    }

    /**
     * Compute the cube surface area.
     *
     * <p>Surface area is useful for shading approximations, texture planning,
     * and physics calculations like drag approximations.</p>
     */
    public double surfaceArea() {
        double a = 6.0 * sideLength * sideLength;
        log.info(() -> "Cube3D surface area computed: " + a);
        return a;
    }

    /**
     * Compute the cube total edge length (sometimes called "perimeter" of the wireframe).
     *
     * <p>A cube has 12 edges. This method returns the sum of all edge lengths.
     * In a rigid cube, this is always {@code 12 * sideLength}, even after rotation/translation.</p>
     */
    public double totalEdgeLength() {
        double p = 12.0 * sideLength;
        log.info(() -> "Cube3D total edge length computed: " + p);
        return p;
    }

    /**
     * Return the length of a single edge.
     *
     * <p>This method exists because callers sometimes need edge size without storing a separate
     * "sideLength" reference. It also helps readability in client code.</p>
     */
    public double edgeLength() {
        log.info(() -> "Cube3D edge length returned: " + sideLength);
        return sideLength;
    }

    /**
     * Update the cube size while keeping the same center.
     *
     * <p>This is useful in editors where the user drags a handle to resize the cube.
     * The cube is rebuilt from the new side length and current center.</p>
     *
     * @param newSideLength the new edge length; must be positive
     */
    public void setSideLength(double newSideLength) {
        log.info(() -> "Updating Cube3D sideLength from " + sideLength + " to " + newSideLength);

        if (newSideLength <= 0.0 || !isFinite(newSideLength)) {
            log.severe("setSideLength failed: invalid newSideLength=" + newSideLength);
            throw new IllegalArgumentException("newSideLength must be a finite number > 0");
        }

        this.sideLength = newSideLength;
        rebuildVerticesFromCenter();
        invalidateCache();
    }

    /**
     * Compute the axis-aligned bounding box (AABB) of this cube in world coordinates.
     *
     * <p>Even when the cube is rotated, an AABB is often used by broad-phase collision and
     * spatial partitioning data structures (BVH, uniform grids) because it is cheap to compare.</p>
     *
     * @return a two-element array: [minCorner, maxCorner]
     */
    public Point3D[] computeAabb() {
        log.info("Computing Cube3D AABB.");

        double minX = Double.POSITIVE_INFINITY, minY = Double.POSITIVE_INFINITY, minZ = Double.POSITIVE_INFINITY;
        double maxX = Double.NEGATIVE_INFINITY, maxY = Double.NEGATIVE_INFINITY, maxZ = Double.NEGATIVE_INFINITY;

        for (Point3D v : vertices) {
            minX = Math.min(minX, v.getX());
            minY = Math.min(minY, v.getY());
            minZ = Math.min(minZ, v.getZ());
            maxX = Math.max(maxX, v.getX());
            maxY = Math.max(maxY, v.getY());
            maxZ = Math.max(maxZ, v.getZ());
        }

        Point3D min = new Point3D(minX, minY, minZ);
        Point3D max = new Point3D(maxX, maxY, maxZ);

        log.info(() -> "AABB computed. min=" + min + " max=" + max);
        return new Point3D[]{min, max};
    }

    /**
     * Determine whether a point lies inside (or on the surface of) this cube.
     *
     * <p>This method is implemented using the cube's local basis assumption only if the cube is axis-aligned.
     * Since this Cube3D supports rotation, the robust "point inside rotated cube" check requires a local
     * coordinate frame (orientation). This class stores vertices only, so here we provide a conservative check:
     * point inside the cube's AABB.</p>
     *
     * <p>That makes this method useful as a fast broad-phase test. For a precise narrow-phase test on a rotated
     * cube, you'd typically store an orientation matrix (or quaternion) and transform the point into local space.</p>
     *
     * @param p the point to test
     * @return true if p is inside the cube's AABB
     */
    public boolean containsPointAabb(Point3D p) {
        Objects.requireNonNull(p, "p must not be null");
        Point3D[] aabb = computeAabb();
        Point3D min = aabb[0];
        Point3D max = aabb[1];

        boolean inside = p.getX() >= min.getX() && p.getX() <= max.getX()
                && p.getY() >= min.getY() && p.getY() <= max.getY()
                && p.getZ() >= min.getZ() && p.getZ() <= max.getZ();

        log.info(() -> "containsPointAabb p=" + p + " -> " + inside);
        return inside;
    }

    // -----------------------------
    // Internal helpers
    // -----------------------------

    /**
     * Rebuild vertex positions from the current center and side length.
     *
     * <p>This method defines the cube in a canonical local layout, then places it at {@code center}.
     * Vertex ordering is stable, which makes edge construction deterministic.</p>
     */
    private void rebuildVerticesFromCenter() {
        double h = sideLength / 2.0;

        // Signs: (-,-,-), (+,-,-), (-,+,-), (+,+,-), (-,-,+), (+,-,+), (-,+,+), (+,+,+)
        vertices[0] = new Point3D(center.getX() - h, center.getY() - h, center.getZ() - h);
        vertices[1] = new Point3D(center.getX() + h, center.getY() - h, center.getZ() - h);
        vertices[2] = new Point3D(center.getX() - h, center.getY() + h, center.getZ() - h);
        vertices[3] = new Point3D(center.getX() + h, center.getY() + h, center.getZ() - h);

        vertices[4] = new Point3D(center.getX() - h, center.getY() - h, center.getZ() + h);
        vertices[5] = new Point3D(center.getX() + h, center.getY() - h, center.getZ() + h);
        vertices[6] = new Point3D(center.getX() - h, center.getY() + h, center.getZ() + h);
        vertices[7] = new Point3D(center.getX() + h, center.getY() + h, center.getZ() + h);

        log.info("Vertices rebuilt from center/sideLength.");
        invalidateCache();
    }

    /**
     * Build the 12 cube edges from the current vertex list.
     *
     * <p>A cube wireframe is defined by edges between vertices that differ by exactly one sign bit
     * in their canonical ordering. This method encodes those 12 connections as index pairs.</p>
     */
    private List<Line3D> buildEdges() {
        // Each pair (a,b) is an edge between vertices[a] and vertices[b]
        int[][] pairs = new int[][]{
                // bottom face (z -h)
                {0, 1}, {1, 3}, {3, 2}, {2, 0},
                // top face (z +h)
                {4, 5}, {5, 7}, {7, 6}, {6, 4},
                // vertical edges
                {0, 4}, {1, 5}, {2, 6}, {3, 7}
        };

        List<Line3D> edges = new ArrayList<>(12);
        for (int[] pair : pairs) {
            edges.add(new Line3D(vertices[pair[0]], vertices[pair[1]]));
        }

        log.info("Edges built: " + edges.size());
        return edges;
    }

    /**
     * Invalidate cached derived data after a transform.
     *
     * <p>This follows a simple cache-invalidation strategy: keep the code correct first,
     * then optimize if needed. Derived values (like edges) can be rebuilt on demand.</p>
     */
    private void invalidateCache() {
        cachedEdges = null;
    }

    /**
     * Create a 3x3 rotation matrix for an axis rotation.
     *
     * <p>This method is intentionally small and focused. Keeping matrix math isolated makes
     * rotation code easier to test and less error-prone.</p>
     */
    private static double[][] rotationMatrix(Axis axis, double radians) {
        double c = Math.cos(radians);
        double s = Math.sin(radians);

        switch (axis) {
            case X:
                return new double[][]{
                        {1, 0, 0},
                        {0, c, -s},
                        {0, s, c}
                };
            case Y:
                return new double[][]{
                        {c, 0, s},
                        {0, 1, 0},
                        {-s, 0, c}
                };
            case Z:
                return new double[][]{
                        {c, -s, 0},
                        {s, c, 0},
                        {0, 0, 1}
                };
            default:
                // Defensive: should never happen with an enum.
                Logger.getLogger(Cube3D.class.getName()).log(Level.SEVERE, "Unknown axis: {0}", axis);
                throw new IllegalArgumentException("Unknown axis: " + axis);
        }
    }

    /**
     * Rotate a point around a pivot using a precomputed rotation matrix.
     *
     * <p>The algorithm:
     * <ol>
     *   <li>Translate point into pivot-relative coordinates.</li>
     *   <li>Multiply by the rotation matrix.</li>
     *   <li>Translate back to world coordinates.</li>
     * </ol>
     *
     * <p>This is the core inner loop of cube rotation and is O(1) per vertex.</p>
     */
    private static Point3D rotatePoint(Point3D p, Point3D pivot, double[][] R) {
        double x = p.getX() - pivot.getX();
        double y = p.getY() - pivot.getY();
        double z = p.getZ() - pivot.getZ();

        double rx = R[0][0] * x + R[0][1] * y + R[0][2] * z;
        double ry = R[1][0] * x + R[1][1] * y + R[1][2] * z;
        double rz = R[2][0] * x + R[2][1] * y + R[2][2] * z;

        return new Point3D(rx + pivot.getX(), ry + pivot.getY(), rz + pivot.getZ());
    }

    /**
     * Utility: finite-number check to guard against NaN/Infinity propagating into geometry.
     *
     * <p>When geometry goes non-finite, renderers and physics engines can fail in hard-to-debug ways.
     * This method provides consistent validation at the boundary of public APIs.</p>
     */
    private static boolean isFinite(double v) {
        return !Double.isNaN(v) && !Double.isInfinite(v);
    }

    @Override
    public String toString() {
        return "Cube3D{center=" + center + ", sideLength=" + sideLength + "}";
    }

    // ----------------------------------------------------------------------
    // Design Patterns + DSA foundations (documentation)
    // ----------------------------------------------------------------------

    /**
     * Patterns used (and why they matter for DSA foundations):
     *
     * <p><b>1) Composition / Aggregation</b><br>
     * The cube is composed of {@code Point3D} vertices and {@code Line3D} edges. This is a foundational OO
     * technique: build complex structures from smaller data types. In data structures terms, the cube’s
     * wireframe is essentially a small graph: vertices (nodes) + edges (connections).</p>
     *
     * <p><b>2) Static Factory Method</b> ({@link #fromMinMax(Point3D, Point3D)})<br>
     * This pattern provides a named constructor that validates input and clarifies intent. Algorithmically,
     * it enforces invariants early (equal extents, positive size) so later operations can assume correctness.</p>
     *
     * <p><b>3) Cache + Invalidation Strategy</b><br>
     * {@link #getEdges()} lazily constructs edges and caches them until a transform occurs. This mirrors
     * a classic time/space tradeoff in algorithms: pay O(1) for repeated reads after an O(n) build.
     * Invalidation keeps the cache correct after updates.</p>
     *
     * <p><b>4) Encapsulation + Defensive Copies</b><br>
     * {@link #getVertices()} returns an immutable snapshot to protect internal arrays. This is the same
     * principle you apply in data structures: control mutation boundaries so invariants remain true.</p>
     *
     * <p><b>How this demonstrates DSA principles</b><br>
     * Most operations here are linear in the number of vertices (n=8): rotation and translation are O(n),
     * where each vertex update is constant time. Edge building is O(1) with a small constant factor (12 edges).
     * The class also models a typical graphics approach: store a compact representation (vertices in an array),
     * derive secondary structures (edges) when needed, and keep transforms as tight loops over arrays for speed.</p>
     */
    @SuppressWarnings("unused")
    private void designNotes() {
        // Documentation-only method.
    }
}
