#include <string>
#include <cstring>
#include <vector>
#include <cstdio>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <omp.h>

template <typename T>
    using Vector3 = Eigen::Matrix<T, 3, 1>;

// the plane is represented by (x - _p) /dot _normal = 0
template <typename T>
class Plane {
public:
    Plane(Vector3<T> p, Vector3<T> normal) {
        _p = p;
        _normal = normal;
        _normal.normalize();
    }

    Vector3<T>& p() { return _p; }
    Vector3<T>& normal() { return _normal; }

    // return if the point is on plane
    // also fill parameter dist as the signed distance from point to plane
    bool onPlane(Vector3<T> point, T& dist) {
        dist = (point - _p).dot(_normal);
        if (std::fabs(dist) < 1e-6) {
            return true;
        } else {
            return false;
        }
    }

private:
    Vector3<T> _p;
    Vector3<T> _normal;
};

template <typename T>
class Triangle {
public:
    Triangle(Vector3<T> v0, Vector3<T> v1, Vector3<T> v2) {
        _vertices[0] = v0;
        _vertices[1] = v1;
        _vertices[2] = v2;
    }

    Vector3<T>* vertices() { return _vertices; }
    Vector3<T>& vertices(int idx) { return _vertices[idx]; }


    // Implement the function to do intersection between triangle and plane p
    // Input: plane p
    // Output: return pairs for intersections with three edges
    //      - enumerate three edges of the triangle and do intersection individually
    //      - consider the case that no intersection
    //      - consider how to avoid repeated intersection points in returned list
    std::vector<Vector3<T>> IntersectPlane(Plane<T> p) {

        std::vector<Vector3<T>> intersections;
        intersections.clear();

                    int i = 0;
                    int ints_ct = 0;
                    while (ints_ct < 2 && i<=2) {
                            int i1 = i;
                            int i2 = (i + 1)%3; 
                            i++;
                            Vector3<T> p1 = vertices(i1);
                            Vector3<T> p2 = vertices(i2);
                            if (std::fabs((p2 - p1).dot(p.normal())) < 1e-6) {
                                    if (std::fabs((p1 - p.p()).dot(p.normal())) < 1e-6) {
                                            intersections.push_back(p1);
                                            intersections.push_back(p2);
                                            ints_ct += 2;
                                    }
                            }
                            else {
                                    T r1 = (p.p() - p1).dot(p.normal()) / ((p2 - p1).dot(p.normal()));
                                    if (r1 >= 0 && r1 <= 1) {
                                            intersections.push_back(p1 + r1 * (p2 - p1));
                                            ints_ct++;
                                    }
                            }
                    }

                    if (intersections.size() == 2) {
                            Vector3<T> p12 = intersections[0] - intersections[1];
                            if (std::sqrt(p12[0] * p12[0] + p12[1] * p12[1] + p12[2] * p12[2]) < 1e-6) {
                                    intersections.clear();
                            }
                    }

        return intersections;
    }

    // Implement the function to do intersection between triangle and a ray
    // Input: a ray, the ray is represented by an origin position and a direction vector
    // Output: return a real number t, the intersection is origin + dir * t, t = -1 means no intersection
    const T IntersectRay(const Vector3<T>& origin, const Vector3<T>& dir) const{

                    T t = -1.0;

                    Vector3<T> u = _vertices[1]- _vertices[0];  //triangle edge vector u=V1-V0
                    Vector3<T> v= _vertices[2] - _vertices[0];   //triangle edge vector v=V2-V0
                    T nx = u[1]*v[2]-u[2]*v[1];
                    T ny = -(u[0]*v[2]-u[2]*v[0]);
                    T nz = u[0]*v[1]-u[1]*v[0];
                    T n_mag = sqrt(nx * nx + ny * ny + nz * nz);
                    Vector3<T> n(nx/n_mag,ny/n_mag,nz/n_mag); //unit normal vector of the triangle

                    //if not parallel to plane triangle lies in
                    T n_dir= (n).dot(dir);

                    if (std::fabs(n_dir) > -(1e-6)) {
                            T r = (n).dot(_vertices[0] - origin) / n_dir;

                            //if plane intersect with Ray
                            if (r >= -(1e-6)) {
                                    Vector3<T> P = origin + r * dir; //intersection point of ray and plane

                                    //determine if intersection point P lies inside or on the triangle
                                    Vector3<T> w = P - _vertices[0]; //vector w=P-V0
                                    //Barycentric Coordinate 
                                    //P=V0+s(V1-V0)+t(V2-V0)=V0+s*u+tt*v
                                    //w=P-V0=s*u+tt*v
                                    T denominator = ((u).dot(v))*((u).dot(v))-((u).dot(u))*((v).dot(v));
                                    T s = ((u.dot(v)) * (v.dot(w)) - (v.dot(v)) * (u.dot(w))) / denominator;
                                    T tt	= ((u.dot(v)) * (w.dot(u)) - (u.dot(u)) * (w.dot(v))) / denominator;
                                    if (s >= -(1e-6) && tt >= -(1e-6) && s + tt <= 1+(1e-6)) {
                                            t = r;
                                    }
                            }
                    }

                    return t;
    }

private:
    Vector3<T> _vertices[3];
};

//Code from: MIT 6.839 Advanced Computer Graphics
//Read in a .stl file and store the triangles and the triangle normals to 
//"triangles", and "normals"
bool ReadSTL(std::string file_name,
    std::vector<std::vector<Vector3<double>>>& triangles, 
    std::vector<Vector3<double>>& normals) {
        
    FILE* fp = std::fopen(file_name.c_str(), "r");

    if (fp == NULL) {
        printf("No STL file found\n");
        return false;
    }

    triangles.clear();
    normals.clear();

    char input[80];
    for (;;) {
        fscanf(fp, "%s", input);
        if (input == std::string("endsolid")) {
            // reach end of file
            break;
        }
        for (;input != std::string("facet");) {
            fscanf(fp, "%s", input);
        }

        std::vector<Vector3<double>> triangle;
        Vector3<double> normal;
        if (std::is_same<double, float>::value) {
            float nx, ny, nz;
            fscanf(fp, "%s %f %f %f\n", input, &nx, &ny, &nz);
            normal[0] = nx; normal[1] = ny; normal[2] = nz;
        }
        else 
            fscanf(fp, "%s %lf %lf %lf\n", input, &normal[0], &normal[1], &normal[2]);

        fscanf(fp, "%s %s", input, input);

        triangle.clear();
        for (int i = 0;i < 3;++i) {
            Vector3<double> p;
            if (std::is_same<double, float>::value) {
                float px, py, pz;
                fscanf(fp, "%s %f %f %f\n", input, &px, &py, &pz);
                p[0] = px; p[1] = py; p[2] = pz;
            }
            else
                fscanf(fp, "%s %lf %lf %lf\n", input, &p[0], &p[1], &p[2]);
            triangle.push_back(p);
        }
        fscanf(fp, "%s %s", input, input);

        triangles.push_back(triangle);
        normals.push_back(normal);
    }

    fclose(fp);
    return true;
}

//normalize the model so that the model is in [-0.5, 0.5]x[-0.5,0.5]x[-0.5, 0.5]
//and zmin should shift to zmin=-0.5
//light source would then have z=zmid=(zmin+zmax)/2
void normalize_model(std::vector<std::vector<Vector3<double>>>& _triangles){
    double xmin=100000;
    double xmax=-100000;
    double ymin=100000;
    double ymax=-100000;
    double zmin=100000;
    double zmax=-100000;
    
    //loop through all triangles
    for (int ii = 0;ii < (signed int)_triangles.size();ii++) {
        Vector3<double> v0 = _triangles[ii][0];//three vertices of a triangle
        Vector3<double> v1 = _triangles[ii][1];
        Vector3<double> v2 = _triangles[ii][2];
        
        if (v0(0) < xmin) { xmin = v0(0); } 
        if (v1(0) < xmin) { xmin = v1(0); } 
        if (v2(0) < xmin) { xmin = v2(0); }
         
        if (v0(0) > xmax) { xmax = v0(0); } 
        if (v1(0) > xmax) { xmax = v1(0); } 
        if (v2(0) > xmax) { xmax = v2(0); }
        
        if (v0(1) < ymin) { ymin = v0(1); } 
        if (v1(1) < ymin) { ymin = v1(1); } 
        if (v2(1) < ymin) { ymin = v2(1); }
         
        if (v0(1) > ymax) { ymax = v0(1); } 
        if (v1(1) > ymax) { ymax = v1(1); } 
        if (v2(1) > ymax) { ymax = v2(1); }
      
        if (v0(2) < zmin) { zmin = v0(2); } 
        if (v1(2) < zmin) { zmin = v1(2); } 
        if (v2(2) < zmin) { zmin = v2(2); }
        
        if (v0(2) > zmax) { zmax = v0(2); } 
        if (v1(2) > zmax) { zmax = v1(2); } 
        if (v2(2) > zmax) { zmax = v2(2); }
    }
    
    double max_range=xmax-xmin;
    if(ymax-ymin>max_range){max_range=ymax-ymin;}
    if(zmax-zmin>max_range){max_range=zmax-zmin;}
    double xmid=0.5*(xmin+xmax);
    double ymid=0.5*(ymin+ymax);
    
    //loop through all triangles
    for (int ii = 0;ii < (signed int)_triangles.size();ii++) {
        _triangles[ii][0](0)=(_triangles[ii][0](0)-xmid)/max_range;
        _triangles[ii][1](0)=(_triangles[ii][1](0)-xmid)/max_range;
        _triangles[ii][2](0)=(_triangles[ii][2](0)-xmid)/max_range;
        
        _triangles[ii][0](1)=(_triangles[ii][0](1)-ymid)/max_range;
        _triangles[ii][1](1)=(_triangles[ii][1](1)-ymid)/max_range;
        _triangles[ii][2](1)=(_triangles[ii][2](1)-ymid)/max_range;
        
        _triangles[ii][0](2)=(_triangles[ii][0](2)-zmin)/max_range-0.5;
        _triangles[ii][1](2)=(_triangles[ii][1](2)-zmin)/max_range-0.5;
        _triangles[ii][2](2)=(_triangles[ii][2](2)-zmin)/max_range-0.5;
    }
    
}

template <typename T>
class Voxelizer {
public:
    Voxelizer(const std::string& stl_file_name, const T dx=0.01) //dx=0.01 corresponds to about 100*100*100 voxels
        : _dx(dx) {
        // Load triangles from the stl file.
        std::vector<Vector3<T>> normals;
        if (!ReadSTL(stl_file_name, _triangles, normals)) {
                std::cout << "ERROR: cannot read " << stl_file_name << std::endl;
                return;
        }
        normalize_model(_triangles);
        // Compute the bounding box of _triangle and save the results into _pmin.
        _pmin = _triangles[0][0];
        Vector3<T> pmax = _triangles[0][0];
        for (const auto& triangle : _triangles)
                for (const auto& v : triangle) {
                        _pmin = _pmin.cwiseMin(v);
                        pmax = pmax.cwiseMax(v);
                }
        
        model_zmin=_pmin[2];
        model_zmax=pmax[2];
        for (int i = 0; i < 3; ++i) {
                _pmin[i] -= _dx;
                pmax[i] += _dx;
        }
        // Compute the number of voxels along each direction.
        for (int i = 0; i < 3; ++i)
                _nvoxel[i] = static_cast<int>((pmax[i] - _pmin[i]) / _dx) + 1;
        // Initialize the voxel array.
        _voxels = std::vector<std::vector<std::vector<bool>>>(_nvoxel.x(),
                std::vector<std::vector<bool>>(_nvoxel.y(),
                        std::vector<bool>(_nvoxel.z(), false)));
    }

    const Vector3<T> pmin() const { return _pmin; }
    const T dx() const { return _dx; }
    const Vector3<T> pmax() const { return _pmin + Vector3<T>(_nvoxel.x(), _nvoxel.y(), _nvoxel.z()) * _dx; }
    const Vector3<int> voxel_num() const { return _nvoxel; }
    double model_zmin;
    double model_zmax;

    // Fill the _voxels array with the correct flag.
    void AdvancedVoxelization() {
        const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2]; // number of voxels in each dimension

        //shoot rays with origins on the y-z grid, in direction of x. (i.e. loop through y and z)
        Vector3<T> ray_dir(1.0, 0.0, 0.0); //direction of the ray
        T dx_half = 0.5 * dx();


        std::vector<T> **intersection = new std::vector<T>*[ny];
        for (int i = 0;i < ny;i++) {
                intersection[i] = new std::vector<T>[nz];
        }

        //loop through all triangles
        for (int ii = 0;ii < (signed int)_triangles.size();ii++) {
            Vector3<T> v0 = _triangles[ii][0];//three vertices of a triangle
            Vector3<T> v1 = _triangles[ii][1];
            Vector3<T> v2 = _triangles[ii][2];
            Triangle<T> currentTria(v0, v1, v2);
            T vy_min = v0[1]; 
            if (v1(1) < vy_min) { vy_min = v1(1); } 
            if (v2(1) < vy_min) { vy_min = v2(1); }
            T vy_max=v0[1]; 
            if (v1(1) > vy_max) { vy_max = v1(1); } 
            if (v2(1) > vy_max) { vy_max = v2(1); }
            T vz_min=v0[2]; 
            if (v1(2) < vz_min) { vz_min = v1(2); } 
            if (v2(2) < vz_min) { vz_min = v2(2); }
            T vz_max=v0[2];  
            if (v1(2) > vz_max) { vz_max = v1(2); } 
            if (v2(2) > vz_max) { vz_max = v2(2); }


            int j_min = (int)((vy_min-_pmin[1])/dx());
            int k_min = (int)((vz_min -_pmin[2])/dx());
            int j_max = (int)((vy_max - _pmin[1])/dx());
            int k_max = (int)((vz_max - _pmin[2])/dx());


            //for each triangle, loop through the rays that can potentially intersect with it
            for(int j=j_min;j<=j_max;j++){
                T center_y = _pmin[1] + dx() * j + dx_half;
                for(int k=k_min;k<=k_max;k++){
                    T center_z = _pmin[2] + dx() * k + dx_half;
                    Vector3<T> ray_o(_pmin[0], center_y, center_z); //origin of the ray
                    T t = currentTria.IntersectRay(ray_o, ray_dir);
                    if (t != -1) {
                            //store the intersection point's x coordinate for the corresponding ray
                            //intersection[j*ny+k].push_back(_pmin[0] + t);
                            intersection[j][k].push_back(_pmin[0] + t);
                    }
                }
            }

        }


        //Loop through the rays
        #pragma omp parallel for num_threads(16)
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {

                //If intersection exists
                if (intersection[j][k].size() != 0) {
                    //sort intersection points's x coordinate
                    //sort edges on a certain layer based on x, ascending order
                    //use Insertion Sorting algorithm
                    //itemsSorted: Number of items that have been sorted so far.
                    for (int itemsSorted = 1; itemsSorted < (signed int)intersection[j][k].size(); itemsSorted++) {
                        // Assume that items A[0], A[1], ... A[itemsSorted-1] 
                        // have already been sorted.  Insert A[itemsSorted]
                        // into the sorted part of the list.

                        T temp = intersection[j][k][itemsSorted];  // The edge to be inserted.
                        int loc = itemsSorted - 1;  // Start at end of list.

                        while (loc >= 0 && intersection[j][k][loc] > temp) {
                                intersection[j][k][loc + 1] = intersection[j][k][loc]; // Bump item from A[loc] up to loc+1.
                                loc = loc - 1;       // Go on to next location.
                        }

                        intersection[j][k][loc + 1] = temp; // Put temp in last vacated space.
                    }
                }

                //determine inside/outside for all the grids along this ray
                for (int i = 0; i < nx; ++i) {

                    if (intersection[j][k].size() == 0) {
                        _voxels[i][j][k] = false;
                    }
                    else {
                        T center_x = _pmin[0] + dx() * i + dx_half;
                        bool search = true; int ii = 0; int ct = 0;
                        while (search && ii < (signed int)intersection[j][k].size() / 2) {
                            if (center_x >= intersection[j][k][2 * ii] && center_x <= intersection[j][k][2 * ii + 1]) {
                                _voxels[i][j][k] = true;
                                search = false;
                            }
                            ii++;
                        }
                        if (search) {
                            _voxels[i][j][k] = false;
                        }

                    }
                }

            }
        }

    }

    double PI=3.14159265358979323846;

    //Code from: MIT 6.839 Advanced Computer Graphics
    void WriteVoxelToMesh(const std::string& stl_file_name) const {
        const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];
        std::vector<std::vector<Vector3<int>>> faces;
        std::vector<Vector3<int>> corners({
            Vector3<int>(0, 0, 0),
            Vector3<int>(0, 0, 1),
            Vector3<int>(0, 1, 0),
            Vector3<int>(0, 1, 1),
            Vector3<int>(1, 0, 0),
            Vector3<int>(1, 0, 1),
            Vector3<int>(1, 1, 0),
            Vector3<int>(1, 1, 1)
        });
        for (int i = 0; i < nx; ++i)
            for (int j = 0; j < ny; ++j)
                for (int k = 0; k < nz; ++k) {
                    if (!_voxels[i][j][k]) continue;
                    // Check -x direction.
                    Vector3<int> cmin(i, j, k);
                    if (i == 0 || !_voxels[i - 1][j][k]) {
                        faces.push_back({ cmin + corners[0], cmin + corners[1], cmin + corners[3] });
                        faces.push_back({ cmin + corners[0], cmin + corners[3], cmin + corners[2] });
                    }
                    if (i == nx - 1 || !_voxels[i + 1][j][k]) {
                        faces.push_back({ cmin + corners[4], cmin + corners[6], cmin + corners[7] });
                        faces.push_back({ cmin + corners[4], cmin + corners[7], cmin + corners[5] });
                    }
                    if (j == 0 || !_voxels[i][j - 1][k]) {
                        faces.push_back({ cmin + corners[0], cmin + corners[4], cmin + corners[5] });
                        faces.push_back({ cmin + corners[0], cmin + corners[5], cmin + corners[1] });
                    }
                    if (j == ny - 1 || !_voxels[i][j + 1][k]) {
                        faces.push_back({ cmin + corners[2], cmin + corners[3], cmin + corners[7] });
                        faces.push_back({ cmin + corners[2], cmin + corners[7], cmin + corners[6] });
                    }
                    if (k == 0 || !_voxels[i][j][k - 1]) {
                        faces.push_back({ cmin + corners[0], cmin + corners[2], cmin + corners[6] });
                        faces.push_back({ cmin + corners[0], cmin + corners[6], cmin + corners[4] });
                    }
                    if (k == nz - 1 || !_voxels[i][j][k + 1]) {
                        faces.push_back({ cmin + corners[5], cmin + corners[7], cmin + corners[3] });
                        faces.push_back({ cmin + corners[5], cmin + corners[3], cmin + corners[1] });
                    }
                }
        std::ofstream fout(stl_file_name);
        fout << "solid vcg" << std::endl;
        for (const auto& f : faces) {
            std::vector<Vector3<T>> p;
            for (const auto& fi : f) {
                Vector3<T> v = _pmin + fi.cast<T>() * _dx;
                p.push_back(v);
            }
            const Vector3<T> n = (p[1] - p[0]).cross(p[2] - p[1]).normalized();
            fout << "  facet normal " << n.x() << " " << n.y() << " " << n.z() << std::endl;
            fout << "    outer loop" << std::endl;
            for (const auto& v : p) {
                fout << "      vertex " << v.x() << " " << v.y() << " " << v.z() << std::endl;
            }
            fout << "    endloop" << std::endl;
            fout << "  endfacet" << std::endl;
        }
        fout << "endsolid vcg" << std::endl;
    }

    //Code from: MIT 6.839 Advanced Computer Graphics
    void WriteVoxelToFile(const std::string &voxel_file) const {
        // File format: _pmin.x _pmin.y _pmin.z dx nx ny nz
        // Then a [nx+1][ny+1][nz+1] array of 0s and 1s.
        const int nx = _nvoxel.x(), ny = _nvoxel.y(), nz = _nvoxel.z();

        // Write it to file.
        std::ofstream fout(voxel_file);
        fout << _pmin.x() << " " << _pmin.y() << " " << _pmin.z() << " " << _dx
                << " " << nx << " " << ny << " " << nz << std::endl;

        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz; ++k) {
                    if (!_voxels[i][j][k]) {
                        fout << "0";
                    }
                    else {
                        fout << "1";
                    }
                }
                fout << "0" << std::endl;
            }
            fout << std::string(nz + 1, '0') << std::endl;
        }
        for (int j = 0; j < ny + 1; ++j) {
            fout << std::string(nz + 1, '0') << std::endl;
        }
    }
    
    
//private:
    std::vector<std::vector<Vector3<T>>> _triangles;
    T _dx;  // The size of each voxel.
    Vector3<T> _pmin;    // The min and max corner of the bounding box.
    Eigen::Vector3i _nvoxel;   // The number of voxels along each direction.
    std::vector<std::vector<std::vector<bool>>> _voxels;   // True <-> voxel is occupied.
};

bool rayThruVoxel(double x0, double y0, double z0, double x1, double y1, double z1, 
        double vxmin, double vxmax, double vymin, double vymax, double vzmin, double vzmax){
    
    double xs=x0; double xb=x1;
    if(x0>x1){xs=x1; xb=x0;}
    double ys=y0; double yb=y1;
    if(y0>y1){ys=y1; yb=y0;}
    double zs=z0; double zb=z1;
    if(z0>z1){zs=z1; zb=z0;}
    
    //fix ray x to be vxmin, find the point (vxmin, y, z) on ray, and see if it goes through voxel
    if((vxmin>=xs && vxmin<=xb)){
        double y_temp=(y1-y0)/(x1-x0)*(vxmin-x0)+y0;
        double z_temp=(z1-z0)/(x1-x0)*(vxmin-x0)+z0;
        if(y_temp>vymin && y_temp<vymax && z_temp>vzmin && z_temp<vzmax){
            return true;
        }
    }
    //fix ray x to be vxmax
    if((vxmax>=xs && vxmax<=xb)){
        double y_temp=(y1-y0)/(x1-x0)*(vxmax-x0)+y0;
        double z_temp=(z1-z0)/(x1-x0)*(vxmax-x0)+z0;
        if(y_temp>vymin && y_temp<vymax && z_temp>vzmin && z_temp<vzmax){
            return true;
        }
    }
    //fix ray y to be vymin
    if((vymin>=ys && vymin<=yb)){
        double x_temp=(x1-x0)/(y1-y0)*(vymin-y0)+x0;
        double z_temp=(z1-z0)/(x1-x0)*(x_temp-x0)+z0;
        if(x_temp>vxmin && x_temp<vxmax && z_temp>vzmin && z_temp<vzmax){
            return true;
        }
    }
    //fix ray y to be vymax
    if((vymax>=ys && vymax<=yb)){
        double x_temp=(x1-x0)/(y1-y0)*(vymax-y0)+x0;
        double z_temp=(z1-z0)/(x1-x0)*(x_temp-x0)+z0;
        if(x_temp>vxmin && x_temp<vxmax && z_temp>vzmin && z_temp<vzmax){
            return true;
        }
    }
    //fix ray z to be vzmin
    if((vzmin>=zs && vzmin<=zb)){
        double x_temp=(x1-x0)/(z1-z0)*(vzmin-z0)+x0;
        double y_temp=(y1-y0)/(x1-x0)*(x_temp-x0)+y0;
        if(x_temp>vxmin && x_temp<vxmax && y_temp>vymin && y_temp<vymax){
            return true;
        }
    }
    //fix ray z to be vzmax
    if((vzmax>=zs && vzmax<=zb)){
        double x_temp=(x1-x0)/(z1-z0)*(vzmax-z0)+x0;
        double y_temp=(y1-y0)/(x1-x0)*(x_temp-x0)+y0;
        if(x_temp>vxmin && x_temp<vxmax && y_temp>vymin && y_temp<vymax){
            return true;
        }
    }
    
    return false;
}

void read_image(const char *fp, int pixel_m, int pixel_n, int **image_grid){
    int i;
    //pixel_m is pixel number in vertical direction; number of rows
    //pixel_n is pixel number in horizontal direction; number of columns
    *image_grid=new int[pixel_m*pixel_n]; 
    //read in silhouette image, with 1-inside, 2-edge, 0-outside of geometry
    std::ifstream file (fp, std::ifstream::in);
    //std::ifstream file ("C:/Users/jiayi/Desktop/Kay/projects/CR_Voro/luckyCloverBinary720_012.txt", std::ifstream::in);
    for (i=0; i<pixel_m; i++) {
        for (int j=0; j<pixel_n; j++){
            file >> (*image_grid)[pixel_n*i+j];
            if((*image_grid)[pixel_n*i+j]==1){
            }
        }
    }
}

//Code from: MIT 6.839 Advanced Computer Graphics
//Marching Cube algorithm
static const int ref_edge_table[256] = {
    0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
};

static const int ref_tri_table[256][16] = {
    { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1 },
    { 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1 },
    { 3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1 },
    { 3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1 },
    { 9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1 },
    { 9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1 },
    { 2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1 },
    { 8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1 },
    { 9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1 },
    { 4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1 },
    { 3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1 },
    { 4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1 },
    { 4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1 },
    { 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1 },
    { 5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1 },
    { 2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1 },
    { 9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1 },
    { 2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1 },
    { 10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1 },
    { 4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1 },
    { 5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1 },
    { 5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1 },
    { 9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1 },
    { 10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1 },
    { 8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1 },
    { 2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1 },
    { 7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1 },
    { 9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1 },
    { 2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1 },
    { 11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1 },
    { 9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1 },
    { 5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1 },
    { 11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1 },
    { 11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1 },
    { 9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1 },
    { 5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1 },
    { 2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1 },
    { 5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1 },
    { 6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1 },
    { 3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1 },
    { 6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1 },
    { 5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1 },
    { 10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1 },
    { 6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1 },
    { 8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1 },
    { 7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1 },
    { 3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1 },
    { 5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1 },
    { 0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1 },
    { 9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1 },
    { 8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1 },
    { 5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1 },
    { 0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1 },
    { 6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1 },
    { 10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1 },
    { 10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1 },
    { 8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1 },
    { 1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1 },
    { 3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1 },
    { 0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1 },
    { 10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1 },
    { 3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1 },
    { 6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1 },
    { 9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1 },
    { 8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1 },
    { 3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1 },
    { 6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1 },
    { 10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1 },
    { 10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1 },
    { 2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1 },
    { 7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1 },
    { 7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1 },
    { 2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1 },
    { 1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1 },
    { 11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1 },
    { 8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1 },
    { 0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1 },
    { 7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1 },
    { 10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1 },
    { 2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1 },
    { 6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1 },
    { 7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1 },
    { 2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1 },
    { 10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1 },
    { 10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1 },
    { 0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1 },
    { 7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1 },
    { 6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1 },
    { 8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1 },
    { 9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1 },
    { 6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1 },
    { 4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1 },
    { 10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1 },
    { 8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1 },
    { 1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1 },
    { 8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1 },
    { 10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1 },
    { 4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1 },
    { 10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1 },
    { 5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1 },
    { 11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1 },
    { 9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1 },
    { 6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1 },
    { 7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1 },
    { 3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1 },
    { 7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1 },
    { 9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1 },
    { 3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1 },
    { 6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1 },
    { 9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1 },
    { 1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1 },
    { 4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1 },
    { 7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1 },
    { 6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1 },
    { 3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1 },
    { 0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1 },
    { 6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1 },
    { 0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1 },
    { 11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1 },
    { 6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1 },
    { 5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1 },
    { 9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1 },
    { 1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1 },
    { 10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1 },
    { 0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1 },
    { 5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1 },
    { 10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1 },
    { 11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1 },
    { 9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1 },
    { 7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1 },
    { 2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1 },
    { 8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1 },
    { 9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1 },
    { 9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1 },
    { 1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1 },
    { 9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1 },
    { 9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1 },
    { 5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1 },
    { 0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1 },
    { 10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1 },
    { 2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1 },
    { 0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1 },
    { 0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1 },
    { 9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1 },
    { 5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1 },
    { 3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1 },
    { 5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1 },
    { 8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1 },
    { 9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1 },
    { 1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1 },
    { 3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1 },
    { 4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1 },
    { 9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1 },
    { 11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1 },
    { 11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1 },
    { 2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1 },
    { 9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1 },
    { 3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1 },
    { 1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1 },
    { 4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1 },
    { 4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1 },
    { 3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1 },
    { 3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1 },
    { 0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1 },
    { 9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1 },
    { 1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { 0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
};
//Code from: MIT 6.839 Advanced Computer Graphics
static const int corner_to_edge[8][3] = {
    { 0, 3, 8 },
    { 0, 1, 9 },
    { 1, 2, 10 },
    { 2, 3, 11 },
    { 4, 7, 8 },
    { 4, 5, 9 },
    { 5, 6, 10 },
    { 6, 7, 11 }
};
//Code from: MIT 6.839 Advanced Computer Graphics
static const int edge_to_corner[12][2] = {
    { 0, 1 },
    { 1, 2 },
    { 2, 3 },
    { 0, 3 },
    { 4, 5 },
    { 5, 6 },
    { 6, 7 },
    { 4, 7 },
    { 0, 4 },
    { 1, 5 },
    { 2, 6 },
    { 3, 7 }
};
//Code from: MIT 6.839 Advanced Computer Graphics
static const int corner_to_xyz[8][3] = {
    { 0, 0, 0 },
    { 0, 0, 1 },
    { 0, 1, 1 },
    { 0, 1, 0 },
    { 1, 0, 0 },
    { 1, 0, 1 },
    { 1, 1, 1 },
    { 1, 1, 0 },
};
//Code from: MIT 6.839 Advanced Computer Graphics
template<typename T>
class MarchingCube {
public:
    MarchingCube(const std::string& voxel_file_name) {
        std::cout << "Entering marching cube..." << std::endl;
        std::ifstream fin(voxel_file_name);
        fin >> _pmin.x() >> _pmin.y() >> _pmin.z();
        fin >> _dx;
        fin >> _nx >> _ny >> _nz;
        std::cout << "pmin = " << _pmin.transpose() << " dx = " << _dx
            << " nx = " << _nx << " ny = " << _ny << " nz = " << _nz << std::endl;
        _corner_flags.resize(_nx + 1);
        for (int i = 0; i < _nx + 1; ++i) {
            _corner_flags[i].resize(_ny + 1);
            for (int j = 0; j < _ny + 1; ++j) {
                fin >> _corner_flags[i][j];
            }
        }
        _triangles.clear();
        GenerateEdgeTable();
    }


    void BuildMesh() {
        _triangles.clear();
        for (int i = 0; i < _nx; ++i)
            for (int j = 0; j < _ny; ++j)
                for (int k = 0; k < _nz; ++k) {
                    // Collect the corners.
                    int cube_idx = 0;
                    for (int c = 0; c < 8; ++c) {
                        if (_corner_flags[i + corner_to_xyz[c][0]][j + corner_to_xyz[c][1]][k + corner_to_xyz[c][2]] == '1')
                            cube_idx |= (1 << c);
                    }
                    for (int t = 0; t < 16; t += 3) {
                        if (ref_tri_table[cube_idx][t] == -1) break;
                        std::vector<Vector3<T>> p(3);
                        for (int ii = 0; ii < 3; ++ii) {
                            const int vid0 = edge_to_corner[ref_tri_table[cube_idx][t + ii]][0];
                            const int vid1 = edge_to_corner[ref_tri_table[cube_idx][t + ii]][1];
                            const Vector3<T> v0 = _pmin + Vector3<T>(i + corner_to_xyz[vid0][0],
                                j + corner_to_xyz[vid0][1], k + corner_to_xyz[vid0][2]) * _dx;
                            const Vector3<T> v1 = _pmin + Vector3<T>(i + corner_to_xyz[vid1][0],
                                j + corner_to_xyz[vid1][1], k + corner_to_xyz[vid1][2]) * _dx;
                            p[ii] = (v0 + v1) * 0.5;
                        }
                        _triangles.push_back(p);
                    }
                }
    }

    void ExportMeshToFile(const std::string& stl_file_name) const {
        std::ofstream fout(stl_file_name);
        fout << "solid vcg" << std::endl;
        for (const auto& p : _triangles) {
            const Vector3<T> n = (p[1] - p[0]).cross(p[2] - p[1]).normalized();
            fout << "  facet normal " << n.x() << " " << n.y() << " " << n.z() << std::endl;
            fout << "    outer loop" << std::endl;
            for (const auto& v : p) {
                fout << "      vertex " << v.x() << " " << v.y() << " " << v.z() << std::endl;
            }
            fout << "    endloop" << std::endl;
            fout << "  endfacet" << std::endl;
        }
        fout << "endsolid vcg" << std::endl;
    }

private:
    void GenerateEdgeTable() {
        for (int i = 0; i < 256; ++i) {
            int edges = 0;
            for (int e = 0; e < 12; ++e) {
                int v0 = edge_to_corner[e][0], v1 = edge_to_corner[e][1];
                bool b0 = (i & (1 << v0)), b1 = (i & (1 << v1));
                if (b0 != b1) edges |= (1 << e);
            }
            _edge_table[i] = edges;
        }

        // Check correctness.
        for (int i = 0; i < 256; ++i) {
            if (_edge_table[i] != ref_edge_table[i]) {
                std::cout << "ERROR: item " << i << " in the edge table is incorrect: "
                    << _edge_table[i] << " vs. " << ref_edge_table[i] << std::endl;
                exit(0);
            }
        }
        std::cout << "Successfully build the edge table..." << std::endl;
    }

    Vector3<T> _pmin;
    int _nx, _ny, _nz;
    T _dx;
    // _corner_flags[i][j][k] == 1 <=> _pmin + (i,j,k) * dx is inside.
    std::vector<std::vector<std::string>> _corner_flags;
    std::vector<std::vector<Vector3<T>>> _triangles;
    int _edge_table[256];
};













int main() {
    
    //image grids: front, back, right, left, up, down
    int *fgrid; //5:3
    int *bgrid; //5:3
    int *rgrid; //5:3
    int *lgrid; //5:3
    int *ugrid; //5:5
    int *dgrid; //5:5
    
    //1. model
    //read in stl
    //scale the model 
    //voxelize model 
printf("A.\n");
    Voxelizer<double> voxelizer("shape3.stl", 0.002);  //0.002
printf("B.\n");
    voxelizer.AdvancedVoxelization();
    double dx_half=0.5*voxelizer._dx;
printf("C.\n");
    //2. patterns
    //read in projection images binary pixel patterns in .txt (6 surrounding walls)
    int b_pixel_m=600; int b_pixel_n=1000;
    int f_pixel_m=600; int f_pixel_n=1000;
    int r_pixel_m=600; int r_pixel_n=1000;
    int l_pixel_m=600; int l_pixel_n=1000;
    int u_pixel_m=1000; int u_pixel_n=1000;
    int d_pixel_m=2000; int d_pixel_n=2000;
    
    read_image("santa_deer.txt", b_pixel_m, b_pixel_n, &bgrid);
    read_image("snow_flakes.txt", f_pixel_m, f_pixel_n, &fgrid);
    read_image("r_Tria_pattern.txt", r_pixel_m, r_pixel_n, &rgrid);
    read_image("l_Tria_pattern.txt", l_pixel_m, l_pixel_n, &lgrid);
    read_image("u_Tria_pattern.txt", u_pixel_m, u_pixel_n, &ugrid);
    read_image("d_Tria_pattern.txt", d_pixel_m, d_pixel_n, &dgrid);
    
printf("D.\n");
    //define surrounding walls domain
    double wall_xmin=-2.5; //bgrid 
    double wall_xmax=2.5;  //fgrid
    double wall_xrange=wall_xmax-wall_xmin;
    double wall_ymin=-2.5; //lgrid
    double wall_ymax=2.5;  //rgrid
    double wall_yrange=wall_ymax-wall_ymin;
    double wall_zmin=-2.5; //dgrid: shape2
    //double wall_zmin=-0.5; //dgrid: spherical
    double wall_zmax=2.5;  //ugrid
    double wall_zrange=wall_zmax-wall_zmin;
    
    //3. putting together
    //define a light source location
    double lsx=0.0; 
    double lsy=0.0; 
    double lsz=voxelizer.model_zmin+0.5*(voxelizer.model_zmax-voxelizer.model_zmin); //at zmid: lamp shape 2
    //double lsz=voxelizer.model_zmin+1.0/3.0*(voxelizer.model_zmax-voxelizer.model_zmin); //at 1/3 height: spherical
    //voxel grid (lsi, lsj, lsk) that the light source is in.
    int lsi=(int)((lsx-voxelizer._pmin[0])/voxelizer._dx);
    int lsj=(int)((lsy-voxelizer._pmin[1])/voxelizer._dx);
    int lsk=(int)((lsz-voxelizer._pmin[2])/voxelizer._dx);
    
    int nx = voxelizer._nvoxel[0], ny = voxelizer._nvoxel[1], nz = voxelizer._nvoxel[2];
printf("E.\n");


    //loop through all pattern grids
    //1. bgrid: santa_deer
    //backwall: x<lsx, x ii from 0 to lsi
    #pragma omp parallel for num_threads(16)
    for (int i=0; i<b_pixel_m; i++) {
        for (int j=0; j<b_pixel_n; j++){
            double x = wall_xmin;
            //if light ray
            //match image grids to surrounding walls coordinates
            if(bgrid[b_pixel_n*i+j]==1){
                double y=wall_ymin+wall_yrange/(b_pixel_n-1)*j;
                double z=wall_zmax-wall_zrange/(b_pixel_m-1)*i;
                
                for(int ii=0; ii<=lsi;ii++){
                    double common=voxelizer._pmin[0]+voxelizer._dx*ii+dx_half;
                    for(int subii=0; subii<10;subii++){
                        double xx=common+0.1*voxelizer._dx*subii;
                        double yy=(y-lsy)/(x-lsx)*(xx-lsx)+lsy;
                        double zz=(z-lsz)/(x-lsx)*(xx-lsx)+lsz;
                        int iitemp=(int)((xx-voxelizer._pmin[0])/voxelizer._dx);
                        int jjtemp=(int)((yy-voxelizer._pmin[1])/voxelizer._dx);
                        int kktemp=(int)((zz-voxelizer._pmin[2])/voxelizer._dx);
                        if(iitemp<nx && jjtemp<ny && kktemp<nz && iitemp>=0 && jjtemp>=0 && kktemp>=0){
                            voxelizer._voxels[iitemp][jjtemp][kktemp]=false;
                        }
                    }
                }
            }
        }
    }
printf("E.1\n");
    //2. fgrid: snowflakes
    //frontwall: x>lsx, x ii from lsi to nx
    #pragma omp parallel for num_threads(16)
    for (int i=0; i<f_pixel_m; i++) {
        for (int j=0; j<f_pixel_n; j++){
            double x = wall_xmax; 
            //if light ray
            //match image grids to surrounding walls coordinates
            if(fgrid[f_pixel_n*i+j]==1){
                double y=wall_ymin+wall_yrange/(f_pixel_n-1)*j;
                double z=wall_zmax-wall_zrange/(f_pixel_m-1)*i;
                
                for(int ii=lsi; ii<nx;ii++){
                    double common=voxelizer._pmin[0]+voxelizer._dx*ii+dx_half;
                    for(int subii=0; subii<10;subii++){
                        double xx=common+0.1*voxelizer._dx*subii;
                        double yy=(y-lsy)/(x-lsx)*(xx-lsx)+lsy;
                        double zz=(z-lsz)/(x-lsx)*(xx-lsx)+lsz;
                        int iitemp=(int)((xx-voxelizer._pmin[0])/voxelizer._dx);
                        int jjtemp=(int)((yy-voxelizer._pmin[1])/voxelizer._dx);
                        int kktemp=(int)((zz-voxelizer._pmin[2])/voxelizer._dx);
                        if(iitemp<nx && jjtemp<ny && kktemp<nz && iitemp>=0 && jjtemp>=0 && kktemp>=0){
                            voxelizer._voxels[iitemp][jjtemp][kktemp]=false;
                        }
                    }
                }
            }
        }
    }
printf("E.2\n");
    //3. lgrid: tria
    //leftwall: y<lsy, y jj from 0 to lsj
    #pragma omp parallel for num_threads(16)
    for (int i=0; i<l_pixel_m; i++) {
        for (int j=0; j<l_pixel_n; j++){
            double y = wall_ymin;
            //if light ray
            //match image grids to surrounding walls coordinates
            if(lgrid[l_pixel_n*i+j]==1){
                double x=wall_xmax-wall_xrange/(l_pixel_n-1)*j;
                double z=wall_zmax-wall_zrange/(l_pixel_m-1)*i;
                
                for(int jj=0; jj<=lsj;jj++){
                    double common=voxelizer._pmin[1]+voxelizer._dx*jj+dx_half;
                    
                    for(int subjj=0; subjj<10;subjj++){
                        
                        double yy=common+0.1*voxelizer._dx*subjj;
                        double xx=(x-lsx)/(y-lsy)*(yy-lsy)+lsx;
                        double zz=(z-lsz)/(x-lsx)*(xx-lsx)+lsz;
                        int iitemp=(int)((xx-voxelizer._pmin[0])/voxelizer._dx);
                        int jjtemp=(int)((yy-voxelizer._pmin[1])/voxelizer._dx);
                        int kktemp=(int)((zz-voxelizer._pmin[2])/voxelizer._dx);
                        if(iitemp<nx && jjtemp<ny && kktemp<nz && iitemp>=0 && jjtemp>=0 && kktemp>=0){
                            
                            voxelizer._voxels[iitemp][jjtemp][kktemp]=false;
                        }
                    }
                }
            }
        }
    }
printf("E.3\n");
    //4. rgrid: tria
    //rightwall: y>lsy, y jj from lsj to ny
    #pragma omp parallel for num_threads(16)
    for (int i=0; i<r_pixel_m; i++) {
        for (int j=0; j<r_pixel_n; j++){
            double y = wall_ymax;
            //if light ray
            //match image grids to surrounding walls coordinates
            if(rgrid[r_pixel_n*i+j]==1){
                double x=wall_xmin+wall_xrange/(r_pixel_n-1)*j;
                double z=wall_zmax-wall_zrange/(r_pixel_m-1)*i;
                
                for(int jj=lsj; jj<ny;jj++){
                    double common=voxelizer._pmin[1]+voxelizer._dx*jj+dx_half;
                    for(int subjj=0; subjj<10;subjj++){
                        double yy=common+0.1*voxelizer._dx*subjj;
                        double xx=(x-lsx)/(y-lsy)*(yy-lsy)+lsx;
                        double zz=(z-lsz)/(x-lsx)*(xx-lsx)+lsz;
                        int iitemp=(int)((xx-voxelizer._pmin[0])/voxelizer._dx);
                        int jjtemp=(int)((yy-voxelizer._pmin[1])/voxelizer._dx);
                        int kktemp=(int)((zz-voxelizer._pmin[2])/voxelizer._dx);
                        if(iitemp<nx && jjtemp<ny && kktemp<nz && iitemp>=0 && jjtemp>=0 && kktemp>=0){
                            voxelizer._voxels[iitemp][jjtemp][kktemp]=false;
                        }
                    }
                }
            }
        }
    }
printf("E.4\n");
    //5. dgrid: tria
    //downwall: z<lsk, z kk from 0 to lsk
    #pragma omp parallel for num_threads(16)
    for (int i=0; i<d_pixel_m; i++) {
        for (int j=0; j<d_pixel_n; j++){
            double z = wall_zmin; 
            //if light ray
            //match image grids to surrounding walls coordinates
            if(dgrid[d_pixel_n*i+j]==1){
                double x=wall_xmin+wall_xrange/(d_pixel_m-1)*i;
                double y=wall_ymin+wall_yrange/(d_pixel_n-1)*j;
                for(int kk=0; kk<=lsk;kk++){
                    double common=voxelizer._pmin[2]+voxelizer._dx*kk+dx_half;
                    for(int subkk=0; subkk<20;subkk++){
                        double zz=common+0.05*voxelizer._dx*subkk;
                        double xx=(x-lsx)/(z-lsz)*(zz-lsz)+lsx;
                        double yy=(y-lsy)/(x-lsx)*(xx-lsx)+lsy;
                        
                        int iitemp=(int)((xx-voxelizer._pmin[0])/voxelizer._dx);
                        int jjtemp=(int)((yy-voxelizer._pmin[1])/voxelizer._dx);
                        int kktemp=(int)((zz-voxelizer._pmin[2])/voxelizer._dx);

                        if(iitemp<nx && jjtemp<ny && kktemp<nz && iitemp>=0 && jjtemp>=0 && kktemp>=0){
                            voxelizer._voxels[iitemp][jjtemp][kktemp]=false;
                        }
                    }
                }
            }
        }
    }
    //6. ugrid: tria
    //upwall: z>lsk, z kk from lsk to nz
    #pragma omp parallel for num_threads(16)
    for (int i=0; i<u_pixel_m; i++) {
        for (int j=0; j<u_pixel_n; j++){
            double z = wall_zmax; 
            //if light ray
            //match image grids to surrounding walls coordinates
            if(ugrid[u_pixel_n*i+j]==1){
                double x=wall_xmin+wall_xrange/(u_pixel_m-1)*i;
                double y=wall_ymax-wall_yrange/(u_pixel_n-1)*j;
                
                for(int kk=lsk; kk<nz;kk++){
                    double common=voxelizer._pmin[2]+voxelizer._dx*kk+dx_half;
                    for(int subkk=0; subkk<10;subkk++){
                        double zz=common+0.1*voxelizer._dx*subkk;
                        double xx=(x-lsx)/(z-lsz)*(zz-lsz)+lsx;
                        double yy=(y-lsy)/(x-lsx)*(xx-lsx)+lsy;
                        int iitemp=(int)((xx-voxelizer._pmin[0])/voxelizer._dx);
                        int jjtemp=(int)((yy-voxelizer._pmin[1])/voxelizer._dx);
                        int kktemp=(int)((zz-voxelizer._pmin[2])/voxelizer._dx);
                        if(iitemp<nx && jjtemp<ny && kktemp<nz && iitemp>=0 && jjtemp>=0 && kktemp>=0){
                            voxelizer._voxels[iitemp][jjtemp][kktemp]=false;
                        }
                    }
                }
            }
        }
    }
    
 
printf("F.\n");
    //output the new voxel model
    voxelizer.WriteVoxelToFile("shape3_voxel_info.txt");
printf("G.\n");
    //marching cube version model
    MarchingCube<double> mc("shape3_voxel_info.txt");
    mc.BuildMesh();
    mc.ExportMeshToFile("shape3_bf_mc.stl");
printf("H.\n");
    return 0;
}




