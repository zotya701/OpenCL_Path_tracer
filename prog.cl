float2 rand(float2 seed)
{
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

typedef struct{
    float3 kd,ks,emission,F0;
    float n, shininess, glossiness;
    int type;
} Material;

typedef struct{
    float3 P,D;
} Ray;
Ray cons_Ray(float3 p, float3 d){
    Ray ray; ray.P=p; ray.D=d; return ray;
}

typedef struct{
    float t;
    float3 P,N;
    Material mat;
} Hit;
Hit cons_Hit(float t, float3 p, float3 n){
    Hit hit; hit.t=t; hit.P=p; hit.N=n; return hit;
}
Hit init_Hit(){
    return cons_Hit(-1.0f, (float3)(0.0f, 0.0f, 0.0f), (float3)(1.0f, 0.0f, 0.0f));
}

float3 Fresnel(Hit hit, Ray old_ray){
    float cosa=fmax(0.0f, dot(hit.N, old_ray.D));
    return hit.mat.F0 + (1-hit.mat.F0)*pow(1-cosa, 5);
}

typedef struct{
    float3 r1,r2,r3,N;
    Material mat;
} Triangle;

typedef struct{
    float3 eye, lookat, up, right;
    float XM, YM;
} Camera;
Ray camera_get_ray(int id, Camera cam, global float2* rnds){
    rnds[id]=rand(rnds[id]);
    int X=cam.XM;
    int Y=cam.YM;
    float x=id%X+rnds[id].x;
    float y=id/X+rnds[id].y;
    float3 p=cam.lookat + cam.right*(2.0f*x/cam.XM-1) + cam.up*(2*y/cam.YM-1);
    float3 d=normalize(p-cam.eye);
    
    return cons_Ray(cam.eye, d);
}
float3 camera_get_view_dir(Hit hit, Camera cam){
    return normalize((cam.eye-hit.P));
}

Hit triangle_intersect(Triangle tri, Ray ray){
    float3 P=ray.P;
    float3 V=ray.D;
    float3 N=tri.N;
    float t=dot((tri.r1-P),N)/dot(V,N);
    if(t<0){
        return init_Hit();
    }
    float3 p=P+V*t;
    if( dot( cross((tri.r2-tri.r1),(p-tri.r1)) , N) >= 0){
        if( dot( cross((tri.r3-tri.r2),(p-tri.r2)) , N)>=0){
            if( dot( cross((tri.r1-tri.r3),(p-tri.r3)) , N)>=0){
                Hit hit=cons_Hit(t,p,N);
                hit.mat=tri.mat;
                return hit;
            }
        }
    }
    return init_Hit();
}

Hit first_intersect(global const Triangle* tris, const int tris_size, const Ray ray){
    Hit best_hit=init_Hit();
    for(int i=0; i<tris_size; ++i){
        Hit hit=triangle_intersect(tris[i], ray);
        //printf("hits[%03d]=\tt=%06.2f \tP=[%06.2f %06.2f %06.2f] \tN=[%06.2f %06.2f %06.2f]\n\r", i, hit.t, hit.P.x, hit.P.y, hit.P.z, hit.N.x, hit.N.y, hit.N.z);
        if(hit.t>0 && (best_hit.t<0 || hit.t<best_hit.t)){
            best_hit=hit;
        }
    }
    return best_hit;
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

Ray new_ray_diffuse(Hit hit, float2 rnds, Ray old_ray){
    float3 X,Y,Z;
    Y=hit.N;
    //if(dot(old_ray.D,hit.N)>0){
    //    Y=-Y;
    //}
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

float3 light_point(Hit hit, float2 rnds){
    float3 p1=(float3)(300.0f, 999.9f, 300.0f);
    float3 p2=(float3)(700.0f, 999.9f, 300.0f);
    float3 p3=(float3)(300.0f, 999.9f, 700.0f);
    float3 p4=(float3)(700.0f, 999.9f, 700.0f);

    float3 v1=p1+(p2-p1)*rnds.x;
    float3 v2=p3+(p4-p3)*rnds.x;

    return v1+(v2-v1)*rnds.y;
}

void kernel trace_ray(global const Triangle* tris, const int tris_size, global Ray* rays, global float2* rnds, const int iterations, const int current_sample, const Camera cam, global float3* colors){
    int id=get_global_id(1)*get_global_size(0) + get_global_id(0);
    float3 factor_A=(float3)(1.0f, 1.0f, 1.0f);
    float3 factor_B=(float3)(1.0f, 1.0f, 1.0f);
    float3 color=(float3)(0.0f, 0.0f, 0.0f);
    
    if(current_sample==0){
        colors[id]=color;
    }

    for(int current=0; current<iterations; ++current){
        Hit hit=first_intersect(tris, tris_size, rays[id]);

        if(hit.t>0){
            if(dot(rays[id].D,hit.N)>0){                                                                                                // to detect if we approach the triangle fromt he front or from behind
                hit.N=-hit.N;
            }
            if(hit.mat.type==0){
                rnds[id]=rand(rnds[id]);
                float3 light_p=light_point(hit,rnds[id]);                                                                                   // random point on light surface
                Ray shadow_ray=cons_Ray(hit.P+hit.N*0.001f, normalize(light_p-(hit.P+hit.N*0.001f)));                                                 // shadow ray from hit point to the light
                Hit shadow_hit=first_intersect(tris, tris_size, shadow_ray);                                                                // calculate the shadow ray's intersection point
                if(shadow_hit.mat.type==3){                                                                                                 // if it was a light emitter
                    float dist=length(light_p-shadow_ray.P)/500;
                    float att=1.0f/(1.0f + 0.1f*dist + 0.01f*dist*dist);

                    float cos_theta=dot(shadow_ray.D, hit.N);
                    float intensity_diffuse=fmax(0.0f, cos_theta);
                    color=color + shadow_hit.mat.emission*hit.mat.kd*intensity_diffuse*factor_A;

                    float3 halfway=normalize(camera_get_view_dir(hit, cam) + shadow_ray.D);
                    float cos_delta=dot(hit.N, halfway);
                    float intensity_specular=fmax(0.0f, cos_delta);
                    color=color + shadow_hit.mat.emission*hit.mat.ks*pow(intensity_specular, hit.mat.shininess)*factor_B;
                }

                Ray old_ray=rays[id];
                rnds[id]=rand(rnds[id]);
                Ray new_ray=new_ray_diffuse(hit, rnds[id], old_ray);
                rays[id]=new_ray;

                float dist=length(light_p-shadow_ray.P)/500;
                float att=1.0f/(1.0f + 0.1f*dist + 0.01f*dist*dist);

                float cos_theta=dot(new_ray.D, hit.N);
                float intensity_diffuse=fmax(0.0f, cos_theta);
                factor_A=factor_A*(hit.mat.kd*intensity_diffuse)*(5);

                float3 halfway=normalize(camera_get_view_dir(hit, cam) + new_ray.D);
                float cos_delta=dot(hit.N, halfway);
                float intensity_specular=fmax(0.0f, cos_delta);
                factor_B=factor_B*(hit.mat.ks*pow(intensity_specular, hit.mat.shininess))*(5);
            }
            if(hit.mat.type==1){
                Ray old_ray=rays[id];
                rays[id]=new_ray_specular(hit, old_ray);
                factor_A=factor_A*Fresnel(hit, old_ray);
                factor_B=factor_B*Fresnel(hit, old_ray);
            }
            if(hit.mat.type==3){
                Ray old_ray=rays[id];
                rnds[id]=rand(rnds[id]);
                Ray new_ray=new_ray_diffuse(hit, rnds[id], old_ray);
                rays[id]=new_ray;
                color=color + hit.mat.emission*(factor_A + factor_B);
            }
        }
    }
    colors[id]=(colors[id]*current_sample + color)/(current_sample+1);
}

void kernel gen_ray(global Ray* rays, const Camera camera, global float2* rnds){
    int id=get_global_id(0);
    rays[id]=camera_get_ray(id, camera, rnds);
}
