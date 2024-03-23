#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/surface/direction_fields.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include <vector>
#include <unordered_map>
#include <math.h>

using namespace geometrycentral;
using namespace geometrycentral::surface;

// Abbreviation for convenience => hash<T,T'> = std::unordered_map<T,T'>
//                              =>    vctr<T> = std::vector<T>
template<typename Key, typename T>
using hash = std::unordered_map<Key, T>;

template<typename T>
using vctr = std::vector<T>;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// Store number of vertices globally, for convenient access.
int n;

// More globals
hash<Edge,int> Gamma;
hash<Corner,double> c;
hash<Corner,double> u;
hash<Edge,double> omega;
hash<Edge,double> gamma;
hash<Corner,double> v_dot;
hash<Edge,double> s;
vctr<double> sigma;
hash<Corner,double> v;
hash<Corner,double> c_tilde;
hash<Corner,double> w;


// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;

// Input: 1-chain on edges Gamma
// Output: a vector containing each vertex in Gamma not on the boundary of the 1-chain
vctr<Vertex> interiorVertices(hash<Edge, int> Gamma) {
    vctr<Vertex> result;
    hash<Vertex, int> count;
    for (Vertex v : mesh->vertices()) {

    }
    for (Edge e : mesh->edges()) {
        if (Gamma[e] != 0) {
            count[e.firstVertex()] += 1;
            count[e.secondVertex()] += 1;
        }
    }
    for (Vertex v : mesh->vertices()) {
        if (count[v] == 2) {
            result.push_back(v);
        }
    }
    return result;
}


// Input: vertex v and 1-chain Gamma on edges
// Output: halfedge in M with its tail being equal to v and it following the curve Gamma
Halfedge outgoingHalfedgeOnCurve(Vertex v, hash<Edge, int> Gamma) {
    for (Halfedge e : mesh->halfedges()) {
        if (e.tailVertex() == v) {
            Edge ed = e.edge();
            if (Gamma[ed] != 0) {
                return e;
            }
        }
    }
}


// Gets the vertex opposite the given halfedge. So the other vert in the tri
Vertex oppositeVertex(Halfedge ij) {
    return ij.next().tipVertex();
}

//
hash<Corner, double> computeReducedCoordinates(hash<Edge, int> Gamma) {
    hash<Corner, double> c;
    for (Vertex i : interiorVertices(Gamma)) {
        if (!(i.isManifold())) {
            continue;
        }
        Halfedge ij0 = outgoingHalfedgeOnCurve(i, Gamma);
        Halfedge ij = ij0;
        int sum = 0;
        do {
            if (!(ij.edge().isBoundary())) {
                Vertex k = oppositeVertex(ij);
                int jump = ij.orientation() ? Gamma[ij.edge()] : -Gamma[ij.edge()];
                sum += jump;
                // TODO: find right corner, centered at i with j and k as other verts
                c[ij.corner()] = sum;
            }
            ij = ij.next().next().twin();
        }
        while (ij != ij0);
    }
    return c;
}

// TODO: InteriorVertices(𝑀, Γ) — returns the set of vertices which are not interior endpoints of the discrete 1-chain Γ on 𝑀.
vctr<Vertex> endpointsOf(hash<Edge,int> Gamma) {
    vctr<Vertex> result;
    return result;
}

vctr<vctr<double>> buildLaplacian(hash<Edge,int> Gamma) {
   vctr<vctr<double>> L(n, vctr<double>(n, 0));
    for (Face curr_face : mesh->faces()) {
        for (Corner curr_corner : curr_face.adjacentCorners()) {
            Vertex k = curr_corner.vertex();
            Vertex i = curr_corner.halfedge().next().vertex();
            Vertex j = curr_corner.halfedge().next().next().vertex();
            bool skip_corner = false;
            for (Vertex other : endpointsOf(Gamma)) {
                if (other == i || other == j) {
                    skip_corner = true;
                    break;
                }
            }
            if (!skip_corner) {
                int idx_i = i.getIndex();
                int idx_j = j.getIndex();
                double theta = geometry->cornerAngle(curr_corner);
                double temp_laplacian_val = 0.5 * (1 / std::tan(theta));
                L[idx_i][idx_i] += temp_laplacian_val;
                L[idx_j][idx_j] += temp_laplacian_val;
                L[idx_i][idx_j] -= temp_laplacian_val;
                L[idx_j][idx_i] -= temp_laplacian_val;
            }
        }
    }
    return L;
}

// Input: A 1-chain Γ ∈ Z|𝐸 | on a mesh 𝑀 = (𝑉 , 𝐸, 𝐹 ) with corner angles 𝜃 , and reduced coordinates 𝑐 ∈ R|𝐶 | (Section 2.3.1).
// Output: The vector 𝑏 ∈ R|𝑉 ∗ | in Equation 10.
vctr<double> buildJumpLaplaceRHS(hash<Edge,int> Gamma, hash<Corner,double> c) {

    // Initialize the vector b
    vctr<double> b(n, 0);

    // Loop over interior verts of Gamma
    for (Vertex curr_vert : interiorVertices(Gamma)) {

        // Loop over all corners incident to curr_vert
        for (Corner curr_corner : curr_vert.adjacentCorners()) {

            // Get triangle face indices from corner
            Vertex i = curr_corner.vertex();
            Vertex j = curr_corner.halfedge().next().vertex();
            Vertex k = curr_corner.halfedge().next().next().vertex();

            // If j or k in endpoints(Gamma) then skip curr_corner
            bool skip_corner = false;
            for (Vertex endpoint_vert : endpointsOf(Gamma)) {
                if (j == endpoint_vert || k == endpoint_vert) {
                    skip_corner = true;
                    break;
                }
            }
            Corner k_centered_corner = curr_corner.halfedge().next().next().corner();
            Corner j_centered_corner = curr_corner.halfedge().next().corner();

            // Fill out vector b if corner is not to be skipped
            if (!skip_corner) {
                b[i.getIndex()] -= 0.5 * (1 / std::tan(geometry->cornerAngle(k_centered_corner))) * c[curr_corner];
                b[j.getIndex()] += 0.5 * (1 / std::tan(geometry->cornerAngle(k_centered_corner))) * c[curr_corner];
                b[i.getIndex()] -= 0.5 * (1 / std::tan(geometry->cornerAngle(j_centered_corner))) * c[curr_corner];
                b[k.getIndex()] += 0.5 * (1 / std::tan(geometry->cornerAngle(j_centered_corner))) * c[curr_corner];
            }
        }
    }
    return b;
}

// TODO: find a sparse PSD solver to use. Eigen probably has one i can use
double solvePSD(vctr<vctr<double>>, vctr<double>) {
    return 0;
}

// Input: A 1-chain Γ ∈ Z|𝐸 | on a mesh 𝑀 = (𝑉 , 𝐸, 𝐹 ) with corner angles 𝜃 , and reduced coordinates 𝑐 ∈ R|𝐶 | .
// Output: A function 𝑢 ∈ R|𝐶 | defined on corners of 𝑀, where 𝑢 solves Equation 10. Values at corners adjacent to endpoints of Γ are left undefined, to be interpolated using Equation 4.
hash<Corner,double> solveJumpEquation(hash<Edge,int> Gamma, hash<Corner,double> c) {

    // Get the laplacian matrix
    vctr<vctr<double>> L = buildLaplacian(Gamma); 

    // Get the RHS of equation (10)
    vctr<double> b = buildJumpLaplaceRHS(Gamma, c);

    // Solve the linear system
    double u0 = solvePSD(L, b);

    // Initialize u, a hash map on edges to doubles
    hash<Corner,double> u;

    // For every corner c' write u[c'] = u0 + c[c']
    for (Corner curr_corner : mesh->corners()) {
        u[curr_corner] = u0 + c[curr_corner];
    }
    return u;
}

// Input: A 1-chain Γ ∈ Z|𝐸 | , and a function 𝑢 ∈ R|𝐶 | with integer jumps across edges of a mesh 𝑀 = (𝑉 , 𝐸, 𝐹 ).
// Output: The Darboux derivative 𝜔 ∈ R|𝐸 | of 𝑢, as a discrete 1-form on edges of 𝑀 (Section 2.4.2).
hash<Edge,double> darbouxDerivative(hash<Edge,int> Gamma, hash<Corner,double> u) {

    // Initialize omega as floating-point valued hash map on edges
    hash<Edge,double> omega;

    // Loop over all edges in mesh
    for (Edge curr_edge : mesh->edges()) {

        // Loop over endpointsOf(Gamma)
        for (Vertex endpoint_vert : endpointsOf(Gamma)) {

            // If our edge is incident to this set we skip
            if (endpoint_vert == curr_edge.firstVertex() || endpoint_vert == curr_edge.secondVertex()) {
                continue;
            }

            // Otherwise set omega[e] = u_jki - u_ijk
            Vertex k = oppositeVertex(curr_edge.halfedge());
            omega[curr_edge] = u[curr_edge.halfedge().next().corner()] - u[curr_edge.halfedge().corner()];
        }
        
    }
    return omega;
}

// Input: A mesh 𝑀 = (𝑉 , 𝐸, 𝐹 ).
// Output: A sparse matrix 𝑑1 ∈ Z|𝐹 | × |𝐸 | representing the discrete exterior derivative on 1-forms.
vctr<vctr<int>> buildOneFormExteriorDerivative() {

    // Store num edges and faces
    int m = mesh->nEdges();
    int k = mesh->nFaces();

    // Initialize the k x m matrix d_1, the exterior derivative on 1-forms
    vctr<vctr<int>> exterior_derivative(k, vctr<int>(m, 0));

    // Loop over faces
    for (Face curr_face : mesh->faces()) {

        // Loop over each edge in face
        for (Edge curr_edge : curr_face.adjacentEdges()) {

            // Based on the orientation of halfedge we say d1[pqr][ij] = +/- 1
            exterior_derivative[curr_face.getIndex()][curr_edge.getIndex()] = curr_edge.halfedge().orientation() ? 1 : -1;
        }
    }
    return exterior_derivative;
}

// Input: A mesh 𝑀 = (𝑉 , 𝐸, 𝐹 ) with corner angles 𝜃 .
// Output: A sparse diagonal matrix ∗1 ∈ Z|𝐹 | × |𝐸 | representing the Hodge star acting on discrete 1-forms.
vctr<vctr<int>> buildOneFormHodgeStar() {
    
    // Store number of edges in m
    int m = mesh->nEdges();

    // Create an empty m x m int matrix
    vctr<vctr<int>> hodge_star(m, vctr<int>(m, 0));

    // Loop over faces
    for (Face curr_face : mesh->faces()) {

        // Loop over each edge in face
        for (Halfedge curr_halfedge : curr_face.adjacentHalfedges()) {

            // Get angle of corner associated with halfedge
            double theta = geometry->cornerAngle(curr_halfedge.corner());

            // Add value to matrix diagonal
            hodge_star[curr_halfedge.edge().getIndex()][curr_halfedge.edge().getIndex()] += int(0.5 * (1 / std::tan(theta)));
        }
    }
}

// Input: A co-closed 1-form 𝜔 ∈ R|𝐸 | on a mesh 𝑀 = (𝑉 , 𝐸, 𝐹 ) with corner angles 𝜃 .
// Output: A harmonic 1-form 𝛾 ∈ R|𝐸 | .
hash<Edge,double> harmonicComponent(hash<Edge,double> omega) {

    // Create hash map to store the result, a real-valued function on E
    hash<Edge,double> result;

    // Function call to exterior derivative
    vctr<vctr<int>> d_1 = buildOneFormExteriorDerivative();

    // Get the hodge start on 1-forms
    vctr<vctr<int>> star_1 = buildOneFormHodgeStar();

    // Solve the sparse positive semi-definite linear system
    // TODO: fix this
    //double beta_tilde = solvePSD(d_1 star_1^{-1} d_1, d_1 omega);
    //delta_beta = 
    return result;
}

// Input: A 1-chain Γ ∈ Z|𝐸 | on a mesh 𝑀 = (𝑉 , 𝐸, 𝐹 ), residual function 𝑣 ∈ R|𝐶 | , and reduced coordinates 𝑐 ∈ R|𝐶 | associated with Γ.
// Output: Updated reduced coordinates ˜𝑐 encoding new jump constraints for the jump Laplace equation (Section 3.5).
hash<Corner,double> subtractJumpDerivative(hash<Edge,int> Gamma, hash<Corner,double> v, hash<Corner,double> c) {

    hash<Corner,double> updated_reduced_coordinates;

    // Loop over interior verts of Gamma
    for (Vertex interior_vert : interiorVertices(Gamma)) {

        // Skip if vertex is not manifold
        if (!(interior_vert.isManifold())) {
            continue;
        }

        // Loop over incident corners
        for (Corner incident_corner : interior_vert.adjacentCorners()) {
            
            // Loop over endpointsOf(Gamma)
            for (Vertex endpoint_vert : endpointsOf(Gamma)) {

                // If our edge is incident to this set we skip
                if (endpoint_vert == incident_corner.halfedge().edge().firstVertex() || endpoint_vert == incident_corner.halfedge().edge().secondVertex()) {
                    continue;
                }

                // If our edge is a boundary edge we also skip
                if (incident_corner.halfedge().edge().isBoundary()) {
                    continue;
                }

                int temp_index = oppositeVertex(incident_corner.halfedge().twin()).getIndex();

                // TODO, get corners based on indices
                Corner ijk;
                Corner ilj;

                updated_reduced_coordinates[ijk] = c[ijk] - (v[ijk] - v[ilj]);
            }
        }
    }
    return updated_reduced_coordinates;
}

// Input: A harmonic 1-form 𝛾 ∈ R|𝐸 | on a mesh 𝑀 = (𝑉 , 𝐸, 𝐹 ).
// Output: Corner values ˚𝑣 𝑗𝑘𝑖 integrating 𝛾 in each triangle of 𝑀.
hash<Corner,double> integrateLocally(hash<Edge,double> gamma) {

    // Initialize result
    hash<Corner,double> local_integral;

    // Initialize g for internal use
    hash<Halfedge,double> g;

    // Loop over faces
    for (Face curr_face : mesh->faces()) {

        // Get ij and jk halfedges
        Halfedge ij = curr_face.halfedge();
        Halfedge jk = ij.next();

        g[ij] = ij.orientation() ? gamma[ij.edge()] : -gamma[ij.edge()];
        g[jk] = jk.orientation() ? gamma[jk.edge()] : -gamma[jk.edge()];

        local_integral[ij.corner()] = 0;
        local_integral[ij.next().corner()] = g[ij];
        local_integral[ij.next().next().corner()] = g[ij] + g[jk];
    }
    return local_integral;
}

// Input: A value ˚𝑣 𝑗𝑘 𝑖 per corner of a mesh 𝑀 = (𝑉 , 𝐸, 𝐹 ).
// Output: Values 𝑠 ∈ R|𝐸 | that give the jump between locally integrated values across each edge of 𝑀.
hash<Edge,double> computeRelativeJumps(hash<Corner,double> v_dot) {

    // Initialize result
    hash<Edge,double> relative_jumps;
    
    // Loop over boundary edges
    for (Edge curr_edge : mesh->edges()) {

        // Check that its a boundary edge
        if (curr_edge.isBoundary()) {

            // Compute the value on current edge
            relative_jumps[curr_edge] = v_dot[curr_edge.halfedge().corner()] - v_dot[curr_edge.halfedge().next().corner()];
        }
    }

    return relative_jumps;
}

// Input: A value ˚𝑣 𝑗𝑘𝑖 per corner of a mesh 𝑀 = (𝑉 , 𝐸, 𝐹 ), and per-triangle shifts 𝜎 ∈ R|𝐹 | .
// Output: A value 𝑣 𝑗𝑘𝑖 per corner describing the residual function.
hash<Corner,double> recoverSolution(hash<Corner,double> v_dot, vctr<double> sigma) {

    // Initialize result hash map on corners
    hash<Corner,double> solution;

    // Loop over all corners
    for (Corner curr_corner : mesh->corners()) {

        // Get face index
        int temp_index = curr_corner.face().getIndex();

        // Residual function
        solution[curr_corner] = v_dot[curr_corner] + sigma[temp_index];
    }

    return solution;
}

// TODO: get actual 1-chains 
// currently just takes the first 5
// edges of the mesh
hash<Edge, int> initGamma() {
    int edges_to_add = 5;
    hash<Edge, int> result; 
    for (Edge e : mesh->edges()) {
        if (edges_to_add > 0)
            result[e] = 1; 
        else
            break;
        edges_to_add -= 1;
    }
    return result;
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {

  if (ImGui::Button("(1) initialize loop")) {
    Gamma = initGamma();
    std::cout << "Loop initialized!" << std::endl;
    for (Edge e : mesh->edges()) {
        if (Gamma[e] != 0) {
            std::cout << "Gamma[" << e << "]=" << Gamma[e] << std::endl;
        }
    }
  }
  if (ImGui::Button("(2) compute reduced coordinates")) {
    c = computeReducedCoordinates(Gamma);
    std::cout << "Reduced coordinates computed!" << std::endl;
    for (Corner cr : mesh->corners()) {
        if (c[cr] != 0) {
            std::cout << "c[" << cr << "]=" << c[cr] << std::endl;
        }
    }
  }
  if (ImGui::Button("(3) solve jump eq")) {
    u = solveJumpEquation(Gamma, c);
    std::cout << "Initial jump equation solved!" << std::endl;
    for (Corner cr : mesh->corners()) {
        if (u[cr] != 0) {
            std::cout << "u[" << cr << "]=" << u[cr] << std::endl;
        }
    }
  }
  if (ImGui::Button("(4) compute darboux derivative")) {
    omega = darbouxDerivative(Gamma, u);
    std::cout << "Darboux derivative computed!" << std::endl;
    for (Edge e : mesh->edges()) {
        std::cout << "Edge " << e << " has value " << omega[e] << " in function omega : E -> R" << std::endl;
    }
  }
  if (ImGui::Button("(5) get harmonic component")) {
    gamma = harmonicComponent(omega);
  }
  if (ImGui::Button("compute surface winding number fn [under construction]")) {

    // SWN Algorithm!!
    hash<Edge,int> Gamma = initGamma();
    hash<Corner,double> c = computeReducedCoordinates(Gamma);
    hash<Corner,double> u = solveJumpEquation(Gamma, c);
    hash<Edge,double> omega = darbouxDerivative(Gamma, u);
    hash<Edge,double> gamma = harmonicComponent(omega);
    hash<Corner,double> v_dot = integrateLocally(gamma);
    hash<Edge,double> s = computeRelativeJumps(v_dot);
    vctr<double> sigma;
    hash<Corner,double> v = recoverSolution(v_dot, sigma);
    hash<Corner,double> c_tilde = subtractJumpDerivative(Gamma, v, c);
    hash<Corner,double> w = solveJumpEquation(Gamma, c);

    CornerData<double> swn;
    for (Corner curr_corner : mesh->corners()) {
        swn[curr_corner] = w[curr_corner];
    }

    psMesh->addCornerScalarQuantity("winding_num", swn);
  }
}


int main(int argc, char **argv) {

  // Configure the argument parser
  args::ArgumentParser parser("geometry-central & Polyscope example project");
  args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help &h) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Make sure a mesh name was given
  if (!inputFilename) {
    std::cerr << "Please specify a mesh file as argument" << std::endl;
    return EXIT_FAILURE;
  }

  // Initialize polyscope
  polyscope::init();

  // Set the callback function
  polyscope::state::userCallback = myCallback;

  // Load mesh
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(args::get(inputFilename));

  // Register the mesh with polyscope
  psMesh = polyscope::registerSurfaceMesh(
      polyscope::guessNiceNameFromPath(args::get(inputFilename)),
      geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));

  // Set vertex tangent spaces
  geometry->requireVertexTangentBasis();
  VertexData<Vector3> vBasisX(*mesh);
  VertexData<Vector3> vBasisY(*mesh);
  for (Vertex v : mesh->vertices()) {
    vBasisX[v] = geometry->vertexTangentBasis[v][0];
    vBasisY[v] = geometry->vertexTangentBasis[v][1];
  }

  // Set n = |V|
  n = mesh->nVertices();

  auto vField =
      geometrycentral::surface::computeSmoothestVertexDirectionField(*geometry);
  psMesh->addVertexTangentVectorQuantity("VF", vField, vBasisX, vBasisY);

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
