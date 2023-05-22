#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <limits>

using namespace std;

typedef unsigned char RGB[3];

#define EPS 0.00001
#define PI 3.14159265359

struct parser::Vec3f& addVec3f( const struct parser::Vec3f &v1, const struct parser::Vec3f &v2);
struct parser::Vec3f& subVec3f( const struct parser::Vec3f &v1, const struct parser::Vec3f &v2);
struct parser::Vec3f& xProductVec3f(const struct parser::Vec3f &v1, const struct parser::Vec3f &v2);
float dotProductVec3f (const struct parser::Vec3f &v1, const struct parser::Vec3f &v2);
struct parser::Vec3f& sProduct (struct parser::Vec3f v, float f);

void showVec3f (const struct parser::Vec3f & v);
void generateImagePlane (const struct parser::Camera &c, std::vector<struct parser::Vec3f> &imagePlane);
float sphere_ray_intersection (const struct parser::Sphere &sph , const struct parser::Vec3f &s, const struct parser::Vec3f &e , const struct parser::Scene &sc);
float face_ray_intersection (const struct parser::Face &fc, const struct parser::Vec3f &s, const struct parser::Vec3f &e, const struct parser::Scene &sc );

struct parser::Vec3f & nVector(const struct parser::Sphere &sp, const struct parser::Vec3f &p, const struct parser::Scene & sc);
struct parser::Vec3f& nVector(const struct parser::Face &fc, const struct parser::Scene &sc);

bool is_in_shadow(const struct parser::Scene &sc, const struct parser::PointLight pl, const struct parser::Face &fc, const struct parser::Vec3f &pt);
bool is_in_shadow(const struct parser::Scene &sc, const struct parser::PointLight pl, const struct parser::Sphere &sp, const struct parser::Vec3f &pt);

float ptp_dist(const struct parser::Vec3f & p1, const struct parser::Vec3f &p2);
bool are_same(const struct parser::Face &f1 ,const struct parser::Face &f2 );
bool are_same(const struct parser::Sphere &s1, const struct parser::Sphere &s2);
bool are_same(const struct parser::Triangle &t1, const struct parser::Triangle &t2 );
void showFace (const struct parser::Face &fc, const struct parser::Scene &sc);

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.
    //
    // Normally, you would be running your ray tracing
    // code here to produce the desired image.

    const RGB BAR_COLOR[8] =
    {
        { 255, 255, 255 },  // 100% White
        { 255, 255,   0 },  // Yellow
        {   0, 255, 255 },  // Cyan
        {   0, 255,   0 },  // Green
        { 255,   0, 255 },  // Magenta
        { 255,   0,   0 },  // Red
        {   0,   0, 255 },  // Blue
        {   0,   0,   0 },  // Black
    };





    struct parser::Camera cam;

    for (int ccc=0; ccc < scene.cameras.size() ; ccc++){
        cam = scene.cameras[ccc];

        int width = cam.image_width, height = cam.image_height;
        int columnWidth = width / 8;
        unsigned char* image = new unsigned char [width * height * 3];
        std::vector<struct parser::Vec3f> imagePlane;
        generateImagePlane(cam, imagePlane);
        bool isRound=false;
        struct parser::Material closest_face_material;
        struct parser::Face closest_face;
        struct parser::Sphere closest_sphere;
        int image_iterator = 0;
    

        for (int i=0 ; i<height ; i++){
            for (int j=0 ; j<width ; j++){
     
                float t = std::numeric_limits<float>::max();
                float t_temp;
                int index = i*width + j;
                for (int k=0 ; k<scene.spheres.size() ; k++){
                    t_temp = sphere_ray_intersection(scene.spheres[k], imagePlane[index], cam.position, scene);
                    if (t_temp > 0 && t_temp < t) {
                        t = t_temp;
                        closest_sphere=scene.spheres[k];
                        isRound=true;
                        
                    }
                }


                for (int k=0 ; k<scene.meshes.size() ; k++){
                    for (int kk=0 ; kk<scene.meshes[k].faces.size() ; kk++){
                        struct parser::Face fc = scene.meshes[k].faces[kk];
                        t_temp = face_ray_intersection(fc, imagePlane[index], cam.position, scene);
                        if (t_temp > 0 && t_temp < t){
                            t = t_temp;
                            closest_face_material=scene.materials[scene.meshes[k].material_id-1];
                            closest_face=scene.meshes[k].faces[kk];
                            isRound=false;
                        }
                    }
                }


                for (int k=0 ; k<scene.triangles.size() ; k++){
                    struct parser::Face fc = scene.triangles[k].indices;
                    t_temp = face_ray_intersection(fc, imagePlane[index], cam.position, scene);
                    if (t_temp > 0 && t_temp < t){
                        t = t_temp;
                        closest_face_material=scene.materials[scene.triangles[k].material_id-1];
                        closest_face=scene.triangles[k].indices;
                        isRound=false;
                    }
                }

                if (t < 0 || t == std::numeric_limits<float>::max()) {
                    image[image_iterator++] = scene.background_color.x;
                    image[image_iterator++] = scene.background_color.y;
                    image[image_iterator++] = scene.background_color.z;
                }
                else {
                    struct parser::Vec3f s=subVec3f(imagePlane[index],cam.position);
                    struct parser::Vec3f inters_point= addVec3f( cam.position, sProduct(s,t));
                    
                    if(isRound){
                        struct parser::Vec3f ambi,spec,color,diff,light_vec;
                        float csin=0;
                        ambi.x=scene.ambient_light.x*scene.materials[closest_sphere.material_id-1].ambient.x;
                        ambi.y=scene.ambient_light.y*scene.materials[closest_sphere.material_id-1].ambient.y;
                        ambi.z=scene.ambient_light.z*scene.materials[closest_sphere.material_id-1].ambient.z;


                        struct parser::Vec3f normal= nVector(closest_sphere,inters_point,scene);

                        diff.x=0;
                        diff.y=0;
                        diff.z=0;
                        spec.x = spec.y = spec.z =0;
                        for(int lll=0;lll<scene.point_lights.size();lll++){
                            bool in_shadow = is_in_shadow(scene, scene.point_lights[lll], closest_sphere, inters_point);
                            light_vec= subVec3f(inters_point, scene.point_lights[lll].position);
                            csin= -1*dotProductVec3f(light_vec,normal)/sqrt(dotProductVec3f(light_vec,light_vec));
                            if(csin<0)
                                csin=0;
                            if(! in_shadow){
                                diff.x+=scene.materials[closest_sphere.material_id-1].diffuse.x* csin *scene.point_lights[lll].intensity.x/dotProductVec3f(light_vec,light_vec);
                                diff.y+=scene.materials[closest_sphere.material_id-1].diffuse.y* csin *scene.point_lights[lll].intensity.y/dotProductVec3f(light_vec,light_vec);
                                diff.z+=scene.materials[closest_sphere.material_id-1].diffuse.z* csin *scene.point_lights[lll].intensity.z/dotProductVec3f(light_vec,light_vec);

                            }
                        }

                        for(int lll=0;lll<scene.point_lights.size();lll++){
                            bool in_shadow = is_in_shadow(scene, scene.point_lights[lll], closest_sphere, inters_point);
                            light_vec= subVec3f(scene.point_lights[lll].position, inters_point);
                            struct parser::Vec3f w_o =  subVec3f(cam.position, inters_point);
                            struct parser::Vec3f w_plus_w = addVec3f(light_vec, w_o);
                            float norm_hv = sqrt(dotProductVec3f(w_plus_w, w_plus_w));
                            struct parser::Vec3f hh = sProduct(w_plus_w, 1.0 / norm_hv);

                            csin = dotProductVec3f(hh,normal);
                            float phong_exp = scene.materials[closest_sphere.material_id-1].phong_exponent;

                            struct parser::Vec3f light_vec_normal =  sProduct(light_vec, (1.0 / sqrt(dotProductVec3f(light_vec,light_vec))));
                            float theta = acos(dotProductVec3f(normal, light_vec_normal));
                            
                            if ( theta > PI / 2.0) {
                                spec.x = spec.y = spec.z = 0;
                            } else {
                                if (! in_shadow){
                                    spec.x += scene.materials[closest_sphere.material_id-1].specular.x * pow(csin, phong_exp) * scene.point_lights[lll].intensity.x / dotProductVec3f(light_vec,light_vec);
                                    spec.y += scene.materials[closest_sphere.material_id-1].specular.y * pow(csin, phong_exp) * scene.point_lights[lll].intensity.y / dotProductVec3f(light_vec,light_vec);
                                    spec.z += scene.materials[closest_sphere.material_id-1].specular.z * pow(csin, phong_exp) * scene.point_lights[lll].intensity.z / dotProductVec3f(light_vec,light_vec);

                                }
                            }
                        }

                        color.x=diff.x+ambi.x + spec.x;
                        color.y=diff.y+ambi.y + spec.y;
                        color.z=diff.z+ambi.z + spec.z;

                        if(color.x>255)
                            color.x=255;
                        if(color.y>255)
                            color.y=255;
                        if(color.z>255)
                            color.z=255;
                        image[image_iterator++]= floor(color.x);
                        image[image_iterator++] = floor(color.y);
                        image[image_iterator++]=floor(color.z);

                    }else{
                        struct parser::Vec3f ambi,spec,color,diff,light_vec;
                        float csin=0;
                        ambi.x=scene.ambient_light.x*closest_face_material.ambient.x;

                        ambi.y=scene.ambient_light.y*closest_face_material.ambient.y;

                        ambi.z=scene.ambient_light.z*closest_face_material.ambient.z;

                        struct parser::Vec3f normal= nVector(closest_face,scene);
                        diff.x=0;
                        diff.y=0;
                        diff.z=0;
                        spec.x = spec.y = spec.z =0;
                        for(int lll=0;lll<scene.point_lights.size();lll++){
                            bool in_shadow = is_in_shadow(scene, scene.point_lights[lll],closest_face, inters_point);
                            light_vec= subVec3f(inters_point, scene.point_lights[lll].position);
                            csin= -1*dotProductVec3f(light_vec,normal)/sqrt(dotProductVec3f(light_vec,light_vec));
                            if(csin<0)
                                csin=0;

                            if (!in_shadow) {
                                diff.x +=closest_face_material.diffuse.x* csin *scene.point_lights[lll].intensity.x/dotProductVec3f(light_vec,light_vec);
                                diff.y +=closest_face_material.diffuse.y* csin *scene.point_lights[lll].intensity.y/dotProductVec3f(light_vec,light_vec);
                                diff.z +=closest_face_material.diffuse.z* csin *scene.point_lights[lll].intensity.z/dotProductVec3f(light_vec,light_vec);
                            }
                        }

                        for(int lll=0;lll<scene.point_lights.size();lll++){
                            bool in_shadow = is_in_shadow(scene, scene.point_lights[lll],closest_face, inters_point);
                            light_vec= subVec3f(scene.point_lights[lll].position, inters_point);
                            struct parser::Vec3f w_o =  subVec3f(cam.position, inters_point);
                            struct parser::Vec3f w_plus_w = addVec3f(light_vec, w_o);
                            float norm_hv = sqrt(dotProductVec3f(w_plus_w, w_plus_w));
                            struct parser::Vec3f hh = sProduct(w_plus_w, 1.0 / norm_hv);

                            csin = dotProductVec3f(hh,normal);
                            struct parser::Vec3f light_vec_normal =  sProduct(light_vec, (1.0 / sqrt(dotProductVec3f(light_vec,light_vec))));
                            float theta = acos(dotProductVec3f(normal, light_vec_normal));
                            
                            if ( theta > PI / 2.0) {
                                spec.x = spec.y = spec.z = 0;
                            }
                            else {
                                float phong_exp = closest_face_material.phong_exponent;
                                if (! in_shadow){
                                    spec.x += closest_face_material.specular.x * pow(csin, phong_exp) * scene.point_lights[lll].intensity.x / dotProductVec3f(light_vec,light_vec);
                                    spec.y += closest_face_material.specular.y * pow(csin, phong_exp) * scene.point_lights[lll].intensity.y / dotProductVec3f(light_vec,light_vec);
                                    spec.z += closest_face_material.specular.z * pow(csin, phong_exp) * scene.point_lights[lll].intensity.z / dotProductVec3f(light_vec,light_vec);
                                }
                            }
                        }

                        /// Shadowing



                        color.x=diff.x+ambi.x + spec.x;
                        color.y=diff.y+ambi.y + spec.y;
                        color.z=diff.z+ambi.z + spec.z;
                        if(color.x>255)
                            color.x=255;
                        if(color.y>255)
                            color.y=255;
                        if(color.z>255)
                            color.z=255;
                        image[image_iterator++] = round(color.x);
                        image[image_iterator++] = round(color.y);
                        image[image_iterator++] = round(color.z);
                    }        
                }
            }
        }

        write_ppm(cam.image_name.c_str(), image, width, height);
    }

}

struct parser::Vec3f& addVec3f( const struct parser::Vec3f &v1, const struct parser::Vec3f &v2){
    struct parser::Vec3f *res = new struct parser::Vec3f;
    res->x = v1.x + v2.x;
    res->y = v1.y + v2.y;
    res->z = v1.z + v2.z;
    return *res;
}

struct parser::Vec3f& subVec3f( const struct parser::Vec3f &v1, const struct parser::Vec3f &v2){
    struct parser::Vec3f *res = new struct parser::Vec3f;
    res->x = v1.x - v2.x;
    res->y = v1.y - v2.y;
    res->z = v1.z - v2.z;
    return *res;
}

struct parser::Vec3f& xProductVec3f(const struct parser::Vec3f &v1, const struct parser::Vec3f &v2){
    struct parser::Vec3f *res = new struct parser::Vec3f;
    res->x = v1.y*v2.z - v2.y*v1.z;
    res->y = v1.z*v2.x - v1.x*v2.z;
    res->z = v1.x*v2.y - v2.x*v1.y;
    return *res;
}

float dotProductVec3f (const struct parser::Vec3f &v1, const struct parser::Vec3f &v2){
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

struct parser::Vec3f& nVector(const struct parser::Face &fc, const struct parser::Scene &sc)
{
    struct parser::Vec3f *nv = new struct parser::Vec3f;
    struct parser::Vec3f a = sc.vertex_data[fc.v0_id-1];
    struct parser::Vec3f b = sc.vertex_data[fc.v1_id-1];
    struct parser::Vec3f c = sc.vertex_data[fc.v2_id-1];

    struct parser::Vec3f aa = subVec3f(b,a);
    struct parser::Vec3f cc = subVec3f(c,b);
    *nv = xProductVec3f(aa,cc);
    
    *nv = sProduct(*nv, 1.0 / sqrt(dotProductVec3f(*nv, *nv)));

    return *nv;
}

struct parser::Vec3f & nVector(const struct parser::Sphere &sp, const struct parser::Vec3f &p, const struct parser::Scene & sc){
    struct parser::Vec3f *nv = new struct parser::Vec3f;
    struct parser::Vec3f c = sc.vertex_data[sp.center_vertex_id-1];
    *nv = sProduct(subVec3f(p, c), 1.0 / sp.radius);
    return *nv;

}

void showVec3f (const struct parser::Vec3f & v){
    std::cout << "(" << v.x  << ", " << v.y << ", "  << v.z << ")\n";
}

void generateImagePlane (const struct parser::Camera &c, std::vector<struct parser::Vec3f> &imagePlane){
    float px = c.position.x;
    float py = c.position.y;
    float pz = c.position.z;
    float lx = c.near_plane.x;
    float ry = c.near_plane.y;
    float bz = c.near_plane.z;
    float tw = c.near_plane.w;
    struct parser::Vec3f v = c.up;
    struct parser::Vec3f w = sProduct(c.gaze, -1.0);
    struct parser::Vec3f u = xProductVec3f(v, w);

    v = sProduct(v, sqrt(1.0 / dotProductVec3f(v,v)));
    w = sProduct(w, sqrt(1.0 / dotProductVec3f(w,w)));
    u = sProduct(u, sqrt(1.0 / dotProductVec3f(u,u)));

    struct parser::Vec3f imgPlaneTL;
    imgPlaneTL.x = px + c.gaze.x * c.near_distance + v.x * tw + u.x * lx;
    imgPlaneTL.y = py + c.gaze.y * c.near_distance + v.y * tw + u.y * lx;
    imgPlaneTL.z = pz + c.gaze.z * c.near_distance + v.z * tw + u.z * lx;

    struct parser::Vec3f imgPlaneTR;
    imgPlaneTR.x = px + c.gaze.x * c.near_distance + v.x * tw + u.x * ry;
    imgPlaneTR.y = py + c.gaze.y * c.near_distance + v.y * tw + u.y * ry;
    imgPlaneTR.z = pz + c.gaze.z * c.near_distance + v.z * tw + u.z * ry;

    struct parser::Vec3f imgPlaneBL;
    imgPlaneBL.x = px + c.gaze.x * c.near_distance + v.x * bz + u.x * lx;
    imgPlaneBL.y = py + c.gaze.y * c.near_distance + v.y * bz + u.y * lx;
    imgPlaneBL.z = pz + c.gaze.z * c.near_distance + v.z * bz + u.z * lx;

    struct parser::Vec3f imgPlaneBR;
    imgPlaneBR.x = px + c.gaze.x * c.near_distance + v.x * bz + u.x * ry;
    imgPlaneBR.y = py + c.gaze.y * c.near_distance + v.y * bz + u.y * ry;
    imgPlaneBR.z = pz + c.gaze.z * c.near_distance + v.z * bz + u.z * ry;



    struct parser::Vec3f rightUnit, downUnit;
    rightUnit = subVec3f(imgPlaneTR, imgPlaneTL);
    rightUnit = sProduct(rightUnit, 1.0 / c.image_width);
    downUnit = subVec3f(imgPlaneBL, imgPlaneTL);
    downUnit = sProduct(downUnit, 1.0 / c.image_height);



    struct parser::Vec3f startingPoint;
    startingPoint.x = imgPlaneTL.x + rightUnit.x / 2.0 + downUnit.x / 2.0;
    startingPoint.y = imgPlaneTL.y + rightUnit.y / 2.0 + downUnit.y / 2.0;
    startingPoint.z = imgPlaneTL.z + rightUnit.z / 2.0 + downUnit.z / 2.0;    

    for (int i = 0; i<c.image_height ; i++){
        for (int j=0 ; j<c.image_width ; j++){
            struct parser::Vec3f temp = addVec3f(sProduct(downUnit,i), sProduct(rightUnit,j));
            imagePlane.push_back(addVec3f(startingPoint, temp));
        }
    }

}

struct parser::Vec3f& sProduct (struct parser::Vec3f v, float f){

    struct parser::Vec3f *nv = new struct parser::Vec3f;
    nv->x = v.x * f;
    nv->y = v.y * f;
    nv->z = v.z * f;
    return *nv;
}

float sphere_ray_intersection (const struct parser::Sphere &sph , const struct parser::Vec3f &s, const struct parser::Vec3f &e , const struct parser::Scene &sc){
    float t;
    float discriminant, d_dot_e_minus_c, d_dot_d,R;
    struct parser::Vec3f e_minus_c, center, d;
    d = subVec3f(s,e);
    center = sc.vertex_data[sph.center_vertex_id-1];
    R = sph.radius;
    e_minus_c = subVec3f(e,center);
    d_dot_e_minus_c = dotProductVec3f(d, e_minus_c);
    d_dot_d = dotProductVec3f(d,d);
    discriminant = pow((dotProductVec3f(d,e_minus_c)),2.0)-d_dot_d*(dotProductVec3f(e_minus_c, e_minus_c)- R*R);
    if ( discriminant < 0) {
        t = -1.0;
    }
    else {
        t = (-1.0 * dotProductVec3f(d,e_minus_c) - sqrt(discriminant)) / d_dot_d;
    }
    return t;
}

float face_ray_intersection (const struct parser::Face &fc, const struct parser::Vec3f &s, const struct parser::Vec3f &e, const struct parser::Scene &sc ){
    struct parser::Vec3f d = subVec3f(s,e);
    float aa, bb, cc ,dd ,ee ,ff ,gg ,hh ,ii, jj, kk, ll;
    struct parser::Vec3f a, b, c;
    a = sc.vertex_data[fc.v0_id-1];
    b = sc.vertex_data[fc.v1_id-1];
    c = sc.vertex_data[fc.v2_id-1];
    aa = a.x - b.x; bb = a.y - b.y ; cc = a.z - b.z ; dd = a.x - c.x ; ee = a.y -c.y ; ff = a.z - c.z ;
    gg = d.x ; hh = d.y ; ii = d.z;
    jj = a.x - e.x ; kk = a.y - e.y ; ll = a.z - e.z;
    float aakk_jjbb = aa*kk - jj*bb;
    float jjcc_aall = jj*cc - aa*ll;
    float bbll_kkcc = bb*ll - kk*cc;
    float eeii_hhff = ee*ii - hh*ff;
    float ggff_ddii = gg*ff - dd*ii;
    float ddhh_eegg = dd*hh - ee*gg;
    float MM = aa*(eeii_hhff) + bb*(ggff_ddii) + cc*(ddhh_eegg);
    float t = (-1.0) * ( ff*aakk_jjbb + ee*jjcc_aall + dd*bbll_kkcc) / MM;
    if ( t < 0) {
        return -1;
    }

    float gamma = (ii*aakk_jjbb + hh*jjcc_aall + gg*bbll_kkcc) / MM;
    if (gamma < 0 || gamma > 1) {
        return -1;
    }

    float beta = (jj*eeii_hhff + kk*ggff_ddii + ll*ddhh_eegg)/ MM;
    if (beta < 0 || beta > 1-gamma){
        return -1;
    }

    return t;


}

float ptp_dist(const struct parser::Vec3f & p1, const struct parser::Vec3f &p2){

    return sqrt( pow(p1.x-p2.x, 2.0) + pow(p1.y-p2.y,2.0) + pow(p1.z-p2.z, 2.0));
}

bool are_same(const struct parser::Face &f1 ,const struct parser::Face &f2 )
{
    if(f1.v0_id != f2.v0_id) return false;
    if(f1.v1_id != f2.v1_id) return false;
    if(f1.v2_id != f2.v2_id) return false;
    return true;
}

bool are_same(const struct parser::Sphere &s1, const struct parser::Sphere &s2)
{
    if (s1.center_vertex_id != s2.center_vertex_id) return false;
    if (s1.radius !=s2.radius) return false;
    return true;
}

bool are_same(const struct parser::Triangle &t1, const struct parser::Triangle &t2 )
{
    if (t1.indices.v0_id != t2.indices.v0_id) return false;
    if (t1.indices.v1_id != t2.indices.v1_id) return false;
    if (t1.indices.v2_id != t2.indices.v2_id) return false;
    return true;
}

bool is_in_shadow(const struct parser::Scene &sc, const struct parser::PointLight pl, const struct parser::Face &fc, const struct parser::Vec3f &pt)
{
    float eps = sc.shadow_ray_epsilon;
    struct parser::Vec3f normal = nVector(fc, sc);
    struct parser::Vec3f ray_start = addVec3f(pt, sProduct(normal, eps));
    struct parser::Vec3f ray_end = pl.position;
    for (int k=0 ; k < sc.meshes.size() ; k++){
        for (int kk=0 ; kk<sc.meshes[k].faces.size() ; kk++){
            struct parser::Face fc_to_check = sc.meshes[k].faces[kk];
            if (are_same(fc_to_check, fc)) continue; // do not check the face with itself.
            float t = face_ray_intersection(fc_to_check, ray_end, ray_start, sc);
            if ( t>0 && t < std::numeric_limits<float>::max()){ // ray hits an object check if before or after light
                float d1 = ptp_dist(ray_end,ray_start); // distance between object and light;
                struct parser::Vec3f x_n_eps= addVec3f(pt, sProduct(normal, eps));
                struct parser::Vec3f twi = sProduct(subVec3f(ray_end, ray_start), t);
                struct parser::Vec3f st = addVec3f(x_n_eps, twi);
                float d2 = ptp_dist(ray_start, st); //distance between object and intersection one;
                if(d2 < d1) {
                    return true;
                }
            }
        }
    }
    for (int k=0 ; k<sc.triangles.size() ; k++){
        struct parser::Face fc_to_check = sc.triangles[k].indices;
        if (are_same(fc_to_check, fc)) continue; // do not check the face with itself.
        float t = face_ray_intersection(fc_to_check, ray_end, ray_start, sc);
        if ( t>0 && t < std::numeric_limits<float>::max()){ // ray hits an object check if before or after light
            float d1 = ptp_dist(ray_end,ray_start); // distance between object and light;
            struct parser::Vec3f x_n_eps= addVec3f(pt, sProduct(normal, eps));
            struct parser::Vec3f twi = sProduct(subVec3f(ray_end, ray_start), t);
            struct parser::Vec3f st = addVec3f(x_n_eps, twi);
            float d2 = ptp_dist(ray_start, st); //distance between object and intersection one;
            if(d2 < d1) {
                return true;
            }
        }
    }

    for (int k=0 ; k<sc.spheres.size() ; k++){
        struct parser::Sphere sp_to_check = sc.spheres[k];
        float t = sphere_ray_intersection(sp_to_check, ray_end, ray_start, sc);
        if ( t>0 && t < std::numeric_limits<float>::max()){ // ray hits an object check if before or after light
            float d1 = ptp_dist(ray_end,ray_start); // distance between object and light;
            struct parser::Vec3f x_n_eps= addVec3f(pt, sProduct(normal, eps));
            struct parser::Vec3f twi = sProduct(subVec3f(ray_end, ray_start), t);
            struct parser::Vec3f st = addVec3f(x_n_eps, twi);
            float d2 = ptp_dist(ray_start, st); //distance between object and intersection one;
            if(d2 < d1) {
                return true;
            }
        }
    }
    return false;

}

bool is_in_shadow(const struct parser::Scene &sc, const struct parser::PointLight pl, const struct parser::Sphere &sp, const struct parser::Vec3f &pt)
{
    float eps = sc.shadow_ray_epsilon;
    struct parser::Vec3f normal = nVector(sp, pt, sc);
    struct parser::Vec3f ray_start = addVec3f(pt, sProduct(normal, eps));
    struct parser::Vec3f ray_end = pl.position;
    for (int k=0 ; k < sc.meshes.size() ; k++){
        for (int kk=0 ; kk<sc.meshes[k].faces.size() ; kk++){
            struct parser::Face fc_to_check = sc.meshes[k].faces[kk];
            float t = face_ray_intersection(fc_to_check, ray_end, ray_start, sc);
            if ( t>0 && t < std::numeric_limits<float>::max()){ // ray hits an object check if before or after light
                float d1 = ptp_dist(ray_end,ray_start); // distance between object and light;
                struct parser::Vec3f x_n_eps= addVec3f(pt, sProduct(normal, eps));
                struct parser::Vec3f twi = sProduct(subVec3f(ray_end, ray_start), t);
                struct parser::Vec3f st = addVec3f(x_n_eps, twi);
                float d2 = ptp_dist(ray_start, st); //distance between object and intersection one;
                if(d2 < d1) {
                    return true;
                }
            }
        }
    }
    for (int k=0 ; k<sc.triangles.size() ; k++){
        struct parser::Face fc_to_check = sc.triangles[k].indices;
        float t = face_ray_intersection(fc_to_check, ray_end, ray_start, sc);
        if ( t>0 && t < std::numeric_limits<float>::max()){ // ray hits an object check if before or after light
            float d1 = ptp_dist(ray_end,ray_start); // distance between object and light;
            struct parser::Vec3f x_n_eps= addVec3f(pt, sProduct(normal, eps));
            struct parser::Vec3f twi = sProduct(subVec3f(ray_end, ray_start), t);
            struct parser::Vec3f st = addVec3f(x_n_eps, twi);
            float d2 = ptp_dist(ray_start, st); //distance between object and intersection one;
            if(d2 < d1) {
                return true;
            }
        }
    }

    for (int k=0 ; k<sc.spheres.size() ; k++){
        struct parser::Sphere sp_to_check = sc.spheres[k];
        float t = sphere_ray_intersection(sp_to_check, ray_end, ray_start, sc);
        if ( t>0 && t < std::numeric_limits<float>::max()){ // ray hits an object check if before or after light
            float d1 = ptp_dist(ray_end,ray_start); // distance between object and light;
            struct parser::Vec3f x_n_eps= addVec3f(pt, sProduct(normal, eps));
            struct parser::Vec3f twi = sProduct(subVec3f(ray_end, ray_start), t);
            struct parser::Vec3f st = addVec3f(x_n_eps, twi);
            float d2 = ptp_dist(ray_start, st); //distance between object and intersection one;
            if(d2 < d1) {
                return true;
            }
        }
    }
    return false;
}

void showFace (const struct parser::Face &fc, const struct parser::Scene &sc){
    cout << "Face ids:\n";
    cout << "( " << fc.v0_id << ", " << fc.v1_id << ", " << fc.v2_id << " )\n";
    cout << "Vertex coordinates are:\n";
    showVec3f(sc.vertex_data[fc.v0_id]);
    showVec3f(sc.vertex_data[fc.v1_id]);
    showVec3f(sc.vertex_data[fc.v2_id]);
}