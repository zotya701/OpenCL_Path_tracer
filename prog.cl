typedef struct{
    float3 bl,tr;
} BBox;

typedef struct{
    int2 trii;
    float split;
    BBox bbox;
} Node;

typedef struct{
    float3 kd,ks,emission,F0;
    float n, shininess, glossiness;
    int type;
} Material;

typedef struct{
    float3 P,D;
} Ray;

bool BBox_intersection(global const BBox* box, Ray* ray, float* tmin){
    float tx1=(box->bl.x-ray->P.x)*native_recip(ray->D.x);
    float tx2=(box->tr.x-ray->P.x)*native_recip(ray->D.x);

    *tmin=fmin(tx1,tx2);
    float tmax=fmax(tx1,tx2);

    float ty1=(box->bl.y-ray->P.y)*native_recip(ray->D.y);
    float ty2=(box->tr.y-ray->P.y)*native_recip(ray->D.y);

    *tmin=fmax(*tmin, fmin(ty1,ty2));
    tmax=fmin(tmax, fmax(ty1,ty2));

    float tz1=(box->bl.z-ray->P.z)*native_recip(ray->D.z);
    float tz2=(box->tr.z-ray->P.z)*native_recip(ray->D.z);

    *tmin=fmax(*tmin, fmin(tz1, tz2));
    tmax=fmin(tmax, fmax(tz1, tz2));

    return tmax>=*tmin;
}

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

void rand(global float2* seed){
    int const a = 16807; //ie 7**5
    int const m = 2147483647; //ie 2**31-1

    int s=(int)(((*seed).x*2.0f-1.0f)*2147483647.0f);
    s = ((long)(s * a))%m;
    (*seed).x=(s/2147483647.0f+1.0f) / 2.0f;

    s=(int)(((*seed).y*2.0f-1.0f)*2147483647.0f);
    s = ((long)(s * a))%m;
    (*seed).y=(s/2147483647.0f+1.0f) / 2.0f;
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

float3 Fresnel(Hit* hit, Ray* old_ray){
    float cosa=fabs(dot(hit->N, old_ray->D));
    return hit->mat.F0 + ((float3)(1.0f, 1.0f, 1.0f)-hit->mat.F0)*pow(1-cosa, 5);
}

Ray cons_Ray(float3 p, float3 d){
    Ray ray; ray.P=p; ray.D=d; return ray;
}

Hit cons_Hit(float t, float3 p, float3 n, ushort mati){
    Hit hit; hit.t=t; hit.P=p; hit.N=n; hit.mati=mati; return hit;
}

Hit init_Hit(){
    return cons_Hit(-1.0f, (float3)(0.0f, 0.0f, 0.0f), (float3)(0.0f, 0.0f, 0.0f), 0);
}

Ray camera_get_ray(int id, const Camera* cam, float2 rnds){
    int X=cam->XM;
    int Y=cam->YM;
    float x=id%X+rnds.x;
    float y=id/X+rnds.y;
    float3 p=cam->lookat + cam->right*(2.0f*x/cam->XM-1) + cam->up*(2*y/cam->YM-1);
    float3 d=normalize(p-cam->eye);
    
    return cons_Ray(cam->eye, d);
}

float3 camera_view_dir(Hit hit, const Camera cam){
    return normalize((cam.eye-hit.P));
}

Hit triangle_intersect(global const Triangle* tri, const Ray* ray){
    Hit hit=init_Hit();
    float3 P=ray->P;
    float3 V=ray->D;
    float3 N=tri->N;
    float t=dot((tri->r1-P),N)/dot(V,N);
    if(t<0){
        return hit;
    }
    float3 p=P+V*t;
    if( dot( cross((tri->r2-tri->r1),(p-tri->r1)) , N) >= 0){
        if( dot( cross((tri->r3-tri->r2),(p-tri->r2)) , N)>=0){
            if( dot( cross((tri->r1-tri->r3),(p-tri->r3)) , N)>=0){
                return cons_Hit(t,p,N,tri->mati);
            }
        }
    }
    return hit;
}

//Hit first_intersect(global const Triangle* tris, const int tris_size, const Ray ray){
Hit first_intersect(global const Triangle* tris, const int from, const int to, const Ray ray){
    Hit best_hit=init_Hit();
    for(int i=from; i<to; ++i){
        Hit hit=triangle_intersect(&tris[i], &ray);
        if(hit.t>0 && (best_hit.t<0 || hit.t<best_hit.t)){
            best_hit=hit;
        }
    }
    return best_hit;
}

void new_ray_diffuse(Hit* hit, float2 rnds, global Ray* ray){
    float3 X,Y,Z;
    Y=hit->N;
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
    (*ray)=cons_Ray(hit->P+Y*0.001f, new_d);
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
    if(disc>0 && rnd>prob){
        return cons_Ray(hit.P - hit.N*0.001f, normalize(old_ray.D*hit.mat.n + hit.N*(cosa*hit.mat.n - sqrt(disc))));
    }else{
        return cons_Ray(hit.P + hit.N*0.001f, normalize(old_ray.D + hit.N*cosa*2.0f));
    }
}

float4 filmic_tone(float3 c){
    float3 c3=(float3)(fmax(0.0f, c.x-0.004f), fmax(0.0f, c.y-0.004f), fmax(0.0f, c.z-0.004f));
    float3 c2=(c3*(c3*6.2f+0.5f))/(c3*(c3*6.2f+1.7f)+0.06f);
    c=pow(c2, 2.2f);
    return (float4)(c.x, c.y, c.z, 1.0f);
}

void stack_push(int* stack, int* ptr, int val){
    if((*ptr)<300){
        stack[*ptr]=val;
        (*ptr)=(*ptr)+1;
    }
}
int stack_pop(int* stack, int* ptr){
    if((*ptr)>0){
        *ptr=*ptr-1;;
        return stack[*ptr];
    }
    return stack[0];
}

Hit kd_intersect(global const Triangle* tris, global const Node* kd_tree, const int kd_tree_size, const Ray ray ){
    Hit hit=init_Hit();
    Hit best_hit=init_Hit();
    int ptr=1;
    float tmin=999999;;
    float dist=0;;
    int stack[300];
    int stack_ptr=0;
    bool found=false;
    bool fail=false;
    bool empty=false;

    while(!empty && !fail && !found){
        if(ptr<kd_tree_size){
            if(BBox_intersection(&kd_tree[ptr].bbox, &ray, &dist)){
                /*
                if(dist>tmin){
                    if(stack_ptr==0){
                        empty=true;
                    }else{
                        ptr=stack_pop(stack, &stack_ptr);
                    }
                }else 
                */
                if(kd_tree[ptr].trii.x<0){
                    stack_push(stack, &stack_ptr, 2*ptr+1);
                    ptr=2*ptr;
                }else{
                    hit=first_intersect(tris, kd_tree[ptr].trii.x, kd_tree[ptr].trii.y, ray);
                    tmin=hit.t;
                    if(hit.t>0 && (best_hit.t<0 || hit.t<best_hit.t)){
                        best_hit=hit;
                    }
                    if(stack_ptr==0){
                        empty=true;
                    }else{
                        ptr=stack_pop(stack, &stack_ptr);
                    }
                }
            }else if(stack_ptr==0){
                empty=true;
            }else{
                ptr=stack_pop(stack, &stack_ptr);
            }
        }else{
            fail=true;
        }
    }
    return best_hit;
}

void kernel trace_ray(write_only image2d_t tex,
                        global const Triangle* tris,
                        const int tris_size,
                        global const Material* materials,
                        global const Node* kd_tree,
                        const int kd_tree_size,
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
        //Hit hit=first_intersect(tris, 0, tris_size, rays[id]);
        Hit hit=kd_intersect(tris, kd_tree, kd_tree_size, rays[id]);

        if(hit.t>0){
            hit.mat=materials[hit.mati];
            if(dot(rays[id].D,hit.N)>0){                                                                                                            // hence the angle between D and N will always be less than 90 degree
                hit.N=-hit.N;
            }
            if(hit.mat.type==0){                                                                                    // diffuse
                rand(&rnds[id]);
                new_ray_diffuse(&hit, rnds[id], &rays[id]);

                float cos_theta=dot(rays[id].D, hit.N);
                float intensity_diffuse=fmax(0.0f, cos_theta);
                factor_A=factor_A*(hit.mat.kd*intensity_diffuse)*factor_S*factor_R;

                float3 halfway=normalize(camera_view_dir(hit, cam) + rays[id].D);
                float cos_delta=dot(hit.N, halfway);
                float intensity_specular=fmax(0.0f, cos_delta);
                factor_B=factor_B*(hit.mat.ks*pow(intensity_specular, hit.mat.shininess))*factor_S*factor_R;
            }
            if(hit.mat.type==1){                                                                                    // specular
                Ray old_ray=rays[id];
                rays[id]=new_ray_specular(hit, old_ray);
                factor_S=factor_S*Fresnel(&hit, &old_ray);
            }
            if(hit.mat.type==2){                                                                                    // refractive
                Ray old_ray=rays[id];
                rand(&rnds[id]);
                rays[id]=new_ray_refractive(hit, old_ray, in, rnds[id].x);
                factor_R=factor_R*(1-Fresnel(&hit, &old_ray));
                in=!in;
            }
            if(hit.mat.type==3){                                                                                    // emitter
                rand(&rnds[id]);
                new_ray_diffuse(&hit, rnds[id], &rays[id]);
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
                rand(&rnds[id]);
                Ray shadow_ray=cons_Ray(hit.P+hit.N*0.001f, light_dir(hit, lights[i], rnds[id]));
                Hit shadow_hit=first_intersect(tris, tris_size, shadow_ray);
                if(shadow_hit.mat.type==3){                                                                                                         // detect if there aren't any obstacles between the object and the light
                    float cos_theta=dot(shadow_ray.D, hit.N);
                    float intensity_diffuse=fmax(0.0f, cos_theta);
                    color=color + shadow_hit.mat.emission*hit.mat.kd*intensity_diffuse*factor_A*factor_S*factor_R;                                  // diffuse color

                    float3 halfway=normalize(camera_view_dir(hit, &cam) + shadow_ray.D);
                    float cos_delta=dot(hit.N, halfway);
                    float intensity_specular=fmax(0.0f, cos_delta);
                    color=color + shadow_hit.mat.emission*hit.mat.ks*pow(intensity_specular, hit.mat.shininess)*factor_B*factor_S*factor_R;         // specular color
                }
            }
            if(hit.mat.type==0){                                                                                    // diffuse
                Ray old_ray=rays[id];
                rand(&rnds[id]);
                Ray new_ray=new_ray_diffuse(hit, rnds[id], old_ray);
                rays[id]=new_ray;

                float cos_theta=dot(new_ray.D, hit.N);
                float intensity_diffuse=fmax(0.0f, cos_theta);
                //factor_A=factor_A*(hit.mat.kd*intensity_diffuse)*(5);
                factor_A=factor_A*(hit.mat.kd*intensity_diffuse)*(5)*factor_S*factor_R;

                float3 halfway=normalize(camera_view_dir(hit, &cam) + new_ray.D);
                float cos_delta=dot(hit.N, halfway);
                float intensity_specular=fmax(0.0f, cos_delta);
                //factor_B=factor_B*(hit.mat.ks*pow(intensity_specular, hit.mat.shininess))*(5);
                factor_B=factor_B*(hit.mat.ks*pow(intensity_specular, hit.mat.shininess))*(5)*factor_S*factor_R;
            }
            if(hit.mat.type==1){                                                                                    // specular
                Ray old_ray=rays[id];
                rays[id]=new_ray_specular(hit, old_ray);
                factor_S=factor_S*Fresnel(&hit, &old_ray);
            }
            if(hit.mat.type==2){                                                                                    // refractive
                Ray old_ray=rays[id];
                rand(&rnds[id]);
                rays[id]=new_ray_refractive(hit, old_ray, in, rnds[id].x);
                Ray new_ray=rays[id];
                factor_R=factor_R*(1-Fresnel(&hit, &old_ray));
                in=!in;
            }
            if(hit.mat.type==3){                                                                                    // emitter
                Ray old_ray=rays[id];
                rand(&rnds[id]);
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
    rand(&rnds[id]);
    rays[id]=camera_get_ray(id, &camera, rnds[id]);
}
