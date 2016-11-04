typedef struct{
    float3 kd,ks,emission,F0;
    float n, shininess, glossiness;
    int type;
} Material;

typedef struct{
    float3 P,D;
} Ray;

typedef struct{
    float t;
    float3 P,N;
    ushort mati;
    Material mat;
} Hit;

typedef struct{
    float3 r1,r2,r3,N;
    ushort mati;
} Triangle;

typedef struct{
    float3 eye, lookat, up, right;
    float XM, YM;
} Camera;

typedef struct{
    float3 p1,p2,p3,p4;
} Light;

float2 rand(float2 seed){
    int s=(int)((seed.x*2.0f-1.0f)*2147483647.0f);
    int const a = 16807; //ie 7**5
    int const m = 2147483647; //ie 2**31-1
    s = ((long)(s * a))%m;
    seed.x=(s/2147483647.0f+1.0f) / 2.0f;

    s=(int)((seed.y*2.0f-1.0f)*2147483647.0f);
    int const a2 = 16807; //ie 7**5
    int const m2 = 2147483647; //ie 2**31-1
    s = ((long)(s * a2))%m2;
    seed.y=(s/2147483647.0f+1.0f) / 2.0f;

    return seed;
}

void ortho_normal_system(const float3* V1, float3* V2, float3* V3){
    const float E=0.001f;
    float3 v1,v2,v3;
    v1=(*V1);
    v2=(*V2);
    v3=(*V3);
    if(fabs(v1.x)<=E && fabs(v1.z)<=E){
        float length=sqrt(v1.y*v1.y + v1.z*v1.z);
        length=1/length;
        v2.x=0;
        v2.y=-v1.z*length;
        v2.z=v1.y*length;
    }
    else
    {
        float length=sqrt(v1.x*v1.x + v1.z*v1.z);
        length=1/length;
        v2.x=-v1.z*length;
        v2.y=0;
        v2.z=v1.x*length;
    }
    v3=cross(v1,v2);
    (*V2)=v2;
    (*V3)=v3;
}

float3 Fresnel(Hit hit, Ray old_ray){
    float cosa=fabs(dot(hit.N, old_ray.D));
    return hit.mat.F0 + ((float3)(1.0f, 1.0f, 1.0f)-hit.mat.F0)*pow(1-cosa, 5);
}

Ray cons_Ray(float3 p, float3 d){
    Ray ray; ray.P=p; ray.D=d; return ray;
}

Hit cons_Hit(float t, float3 p, float3 n, ushort mati){
    Hit hit; hit.t=t; hit.P=p; hit.N=n; hit.mati=mati; return hit;
}

Hit init_Hit(){
    return cons_Hit(-1.0f, (float3)(0.0f, 0.0f, 0.0f), (float3)(0.0f, 0.0f, 0.0f), -1);
}

Ray camera_get_ray(int id, Camera cam, float2 rnds){
    int X=cam.XM;
    int Y=cam.YM;
    float x=id%X+rnds.x;
    float y=id/X+rnds.y;
    float3 p=cam.lookat + cam.right*(2.0f*x/cam.XM-1) + cam.up*(2*y/cam.YM-1);
    float3 d=normalize(p-cam.eye);
    
    return cons_Ray(cam.eye, d);
}

float3 camera_view_dir(Hit hit, Camera cam){
    return normalize((cam.eye-hit.P));
}

Hit triangle_intersect(Triangle tri, Ray ray){
    Hit hit=init_Hit();
    float3 P=ray.P;
    float3 V=ray.D;
    float3 N=tri.N;
    float t=dot((tri.r1-P),N)/dot(V,N);
    if(t<0){
        return hit;
    }
    float3 p=P+V*t;
    if( dot( cross((tri.r2-tri.r1),(p-tri.r1)) , N) >= 0){
        if( dot( cross((tri.r3-tri.r2),(p-tri.r2)) , N)>=0){
            if( dot( cross((tri.r1-tri.r3),(p-tri.r3)) , N)>=0){
                return cons_Hit(t,p,N,tri.mati);
            }
        }
    }
    return hit;
}

Hit first_intersect(global const Triangle* tris, const int tris_size, const Ray ray){
    Hit best_hit=init_Hit();
    for(int i=0; i<tris_size; ++i){
        Hit hit=triangle_intersect(tris[i], ray);
        if(hit.t>0 && (best_hit.t<0 || hit.t<best_hit.t)){
            best_hit=hit;
        }
    }
    return best_hit;
}

Ray new_ray_diffuse(Hit hit, float2 rnds, Ray old_ray){
    float3 X,Y,Z;
    Y=hit.N;
    ortho_normal_system(&Y,&Z,&X);
    float rnd1,rnd2,r,theta,x,y,z;
    rnd1=rnds.x;
    rnd2=rnds.y;
    r=sqrt(rnd1);
    theta=2*M_PI*rnd2;
    x=r*cos(theta);
    y=r*sin(theta);
    z=sqrt(1-rnd1);
    float3 new_d=normalize(X*x+Y*z+Z*y);
    return cons_Ray(hit.P+Y*0.001f, new_d);
}

Ray new_ray_specular(Hit hit, Ray old_ray){
    float3 new_d=normalize(old_ray.D - hit.N*dot(hit.N, old_ray.D)*2.0f);
    return cons_Ray(hit.P+hit.N*0.001f, new_d);
}

Ray new_ray_refractive(Hit hit, Ray old_ray, bool in, float rnd){
    if(in){
        hit.mat.n=1.0f/hit.mat.n;
    }
    hit.mat.n=1.0f/hit.mat.n;
    float cosa=-dot(old_ray.D, hit.N);
    float disc=1.0f - (1.0f - cosa*cosa)*hit.mat.n*hit.mat.n;
    float3 F=hit.mat.F0 + ((float3)(1.0f, 1.0f, 1.0f) - hit.mat.F0)*pow(1-cosa, 5);
    float prob=(F.x+F.y+F.z)/3.0f;
    if(disc>0 && rnd){
        return cons_Ray(hit.P - hit.N*0.001f, normalize(old_ray.D*hit.mat.n + hit.N*(cosa*hit.mat.n - sqrt(disc))));
    }else{
        return cons_Ray(hit.P + hit.N*0.001f, normalize(old_ray.D + hit.N*cosa*2.0f));
    }
}

float3 light_dir(Hit hit, Light l, float2 rnds){
    float3 v1=l.p1+(l.p2-l.p1)*rnds.x;
    float3 v2=l.p4+(l.p3-l.p4)*rnds.x;
    float3 v=v1+(v2-v1)*rnds.y;
    return normalize(v-(hit.P+hit.N*0.001f));
}

float4 filmic_tone(float3 c){
    float3 c3=(float3)(fmax(0.0f, c.x-0.004f), fmax(0.0f, c.y-0.004f), fmax(0.0f, c.z-0.004f));
    float3 c2=(c3*(c3*6.2f+0.5f))/(c3*(c3*6.2f+1.7f)+0.06f);
    c=pow(c2, 2.2f);
    return (float4)(c.x, c.y, c.z, 1.0f);
}


/*
kd_search( tree, ray ){
    (global_tmin, global_tmax) = intersect( tree.bounds, ray )
    search_node( tree.root, ray, global_tmin, global_tmax )
}
search_node( node, ray, tmin, tmax ){
    if( node.is_leaf ){
        search_leaf( node, ray, tmin, tmax )
    }else{
        search_split( node, ray, tmin, tmax )
    }
}
search_split( split, ray, tmin, tmax ){
    a = split.axis
    thit = ( split.value - ray.origin[a] ) / ray.direction[a]
    (first, second) = order( ray.direction[a], split.left, split.right )
    if( thit >= tmax or thit < 0 ){
        search_node( first, ray, tmin, tmax )
    }else if( thit <= tmin ){
        search_node( second, ray, tmin, tmax )
    }else{
        search_node( first, ray, tmin, thit )
    }
}
search_leaf( leaf, ray, tmin, tmax ){
    // search for a hit in this leaf
    if( found_hit and hit.t < tmax ){
        succeed( hit )
    }else{
        continue_search( leaf, ray, tmin, tmax )
    }
}
continue_search( leaf, ray, tmin, tmax ){
    if( tmax == global_tmax ){
        fail()
    }else{
        tmin = tmax
        tmax = global_tmax
        search_node( tree.root, ray, tmin, tmax )
    }
}
*/


void kernel trace_ray(write_only image2d_t tex,
                        global const Triangle* tris,
                        const int tris_size,
                        global const Material* materials,
                        global const Light* lights,
                        const int lights_size,
                        global Ray* rays,
                        global float2* rnds,
                        const int iterations,
                        const int current_sample,
                        const Camera cam,
                        global float3* colors){
    int id=get_global_id(1)*get_global_size(0) + get_global_id(0);
    float3 factor_A=(float3)(1.0f, 1.0f, 1.0f);
    float3 factor_B=(float3)(1.0f, 1.0f, 1.0f);
    float3 factor_S=(float3)(1.0f, 1.0f, 1.0f);
    float3 factor_R=(float3)(1.0f, 1.0f, 1.0f);
    float3 color=(float3)(0.1f, 0.1f, 0.1f);
    
    if(current_sample==0){
        colors[id]=color;
    }
    
    bool in=false;
    for(int current=0; current<iterations; ++current){
        Hit hit=first_intersect(tris, tris_size, rays[id]);

        if(hit.t>0){
            hit.mat=materials[hit.mati];
            if(dot(rays[id].D,hit.N)>0){                                                                                                            // hence the angle between D and N will always be less than 90 degree
                hit.N=-hit.N;
            }
            if(hit.mat.type==0){                                                                                    // diffuse
                Ray old_ray=rays[id];
                rnds[id]=rand(rnds[id]);
                Ray new_ray=new_ray_diffuse(hit, rnds[id], old_ray);
                rays[id]=new_ray;

                float cos_theta=dot(new_ray.D, hit.N);
                float intensity_diffuse=fmax(0.0f, cos_theta);
                factor_A=factor_A*(hit.mat.kd*intensity_diffuse)*factor_S*factor_R;

                float3 halfway=normalize(camera_view_dir(hit, cam) + new_ray.D);
                float cos_delta=dot(hit.N, halfway);
                float intensity_specular=fmax(0.0f, cos_delta);
                factor_B=factor_B*(hit.mat.ks*pow(intensity_specular, hit.mat.shininess))*factor_S*factor_R;
            }
            if(hit.mat.type==1){                                                                                    // specular
                Ray old_ray=rays[id];
                rays[id]=new_ray_specular(hit, old_ray);
                factor_S=factor_S*Fresnel(hit, old_ray);
            }
            if(hit.mat.type==2){                                                                                    // refractive
                Ray old_ray=rays[id];
                rnds[id]=rand(rnds[id]);
                rays[id]=new_ray_refractive(hit, old_ray, in, rnds[id].x);
                Ray new_ray=rays[id];
                factor_R=factor_R*(1-Fresnel(hit, old_ray));
                in=!in;
            }
            if(hit.mat.type==3){                                                                                    // emitter
                Ray old_ray=rays[id];
                rnds[id]=rand(rnds[id]);
                Ray new_ray=new_ray_diffuse(hit, rnds[id], old_ray);
                rays[id]=new_ray;
                color=color + hit.mat.emission*(factor_A + factor_B)*factor_S*factor_R*(current+1);
            }
        }else{
            break;
        }
    }
    

    /*
    bool in=false;
    for(int current=0; current<iterations; ++current){
        Hit hit=first_intersect(tris, tris_size, rays[id]);

        if(hit.t>0){
            if(dot(rays[id].D,hit.N)>0){                                                                                                            // hence the angle between D and N will always be less than 90 degree
                hit.N=-hit.N;
            }
            for(int i=0;i<lights_size;++i){                                                                                                         // sample all lights
                rnds[id]=rand(rnds[id]);
                Ray shadow_ray=cons_Ray(hit.P+hit.N*0.001f, light_dir(hit, lights[i], rnds[id]));
                Hit shadow_hit=first_intersect(tris, tris_size, shadow_ray);
                if(shadow_hit.mat.type==3){                                                                                                         // detect if there aren't any obstacles between the object and the light
                    float cos_theta=dot(shadow_ray.D, hit.N);
                    float intensity_diffuse=fmax(0.0f, cos_theta);
                    color=color + shadow_hit.mat.emission*hit.mat.kd*intensity_diffuse*factor_A*factor_S*factor_R;                                  // diffuse color

                    float3 halfway=normalize(camera_view_dir(hit, cam) + shadow_ray.D);
                    float cos_delta=dot(hit.N, halfway);
                    float intensity_specular=fmax(0.0f, cos_delta);
                    color=color + shadow_hit.mat.emission*hit.mat.ks*pow(intensity_specular, hit.mat.shininess)*factor_B*factor_S*factor_R;         // specular color
                }
            }
            if(hit.mat.type==0){                                                                                    // diffuse
                Ray old_ray=rays[id];
                rnds[id]=rand(rnds[id]);
                Ray new_ray=new_ray_diffuse(hit, rnds[id], old_ray);
                rays[id]=new_ray;

                float cos_theta=dot(new_ray.D, hit.N);
                float intensity_diffuse=fmax(0.0f, cos_theta);
                //factor_A=factor_A*(hit.mat.kd*intensity_diffuse)*(5);
                factor_A=factor_A*(hit.mat.kd*intensity_diffuse)*(5)*factor_S*factor_R;

                float3 halfway=normalize(camera_view_dir(hit, cam) + new_ray.D);
                float cos_delta=dot(hit.N, halfway);
                float intensity_specular=fmax(0.0f, cos_delta);
                //factor_B=factor_B*(hit.mat.ks*pow(intensity_specular, hit.mat.shininess))*(5);
                factor_B=factor_B*(hit.mat.ks*pow(intensity_specular, hit.mat.shininess))*(5)*factor_S*factor_R;
            }
            if(hit.mat.type==1){                                                                                    // specular
                Ray old_ray=rays[id];
                rays[id]=new_ray_specular(hit, old_ray);
                factor_S=factor_S*Fresnel(hit, old_ray);
            }
            if(hit.mat.type==2){                                                                                    // refractive
                Ray old_ray=rays[id];
                rnds[id]=rand(rnds[id]);
                rays[id]=new_ray_refractive(hit, old_ray, in, rnds[id].x);
                Ray new_ray=rays[id];
                factor_R=factor_R*(1-Fresnel(hit, old_ray));
                in=!in;
            }
            if(hit.mat.type==3){                                                                                    // emitter
                Ray old_ray=rays[id];
                rnds[id]=rand(rnds[id]);
                Ray new_ray=new_ray_diffuse(hit, rnds[id], old_ray);
                rays[id]=new_ray;
                color=color + hit.mat.emission*(factor_A + factor_B)*factor_S*factor_R;
            }
        }
    }
    */
    colors[id]=(colors[id]*current_sample + color)/(current_sample+1);
    write_imagef(tex, (int2)(get_global_id(0), get_global_id(1)), filmic_tone(colors[id]));
}

void kernel gen_ray(global Ray* rays, const Camera camera, global float2* rnds){
    int id=get_global_id(0);
    rnds[id]=rand(rnds[id]);
    rays[id]=camera_get_ray(id, camera, rnds[id]);
}
