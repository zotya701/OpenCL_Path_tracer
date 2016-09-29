
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

typedef struct{
    float3 r1,r2,r3,N;
    Material mat;
} Triangle;

typedef struct{
    float3 eye, lookat, up, right;
    float XM, YM;
} Camera;
Ray camera_get_ray(int id, Camera cam){
    int X=cam.XM;
    int Y=cam.YM;
    float x=id%X+0.5f;
    float y=id/X+0.5f;
    float3 p=cam.lookat + cam.right*(2.0f*x/cam.XM-1) + cam.up*(2*y/cam.YM-1);
    float3 d=normalize(p-cam.eye);
    
    return cons_Ray(cam.eye, d);
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

Ray new_ray_diffuse(Hit hit){
    float3 X,Y,Z;
    Y=hit.N;
    ortho_normal_system(&Y,&Z,&X);
    float rnd1,rnd2,r,theta,x,y,z;
    //rnd1=RND2;
    //rnd2=RND2;
    rnd1=0.0f;
    rnd2=0.0f;
    r=sqrt(rnd1);
    theta=2*M_PI*rnd2;
    x=r*cos(theta);
    y=r*sin(theta);
    z=sqrt(fmax(0.0f, 1-rnd1));
    float3 new_d = normalize(X*x+Y*z+Z*y);
    return cons_Ray(hit.P,new_d);
}

void kernel trace_ray(global const Triangle* tris, const int tris_size, global Ray* rays, global float3* colors){
    int id=get_global_id(0);
    Hit hit=first_intersect(tris, tris_size, rays[id]);

    

    //printf("hits[%03d]=\tt=%06.2f \tP=[%06.2f %06.2f %06.2f] \tN=[%06.2f %06.2f %06.2f]\n\r", id, hit.t, hit.P.x, hit.P.y, hit.P.z, hit.N.x, hit.N.y, hit.N.z);

    Ray old_ray=rays[id];
    Ray new_ray=new_ray_diffuse(hit);
    rays[id]=new_ray;
    if(hit.t>0){
        colors[id]=hit.mat.kd;
    }else{
        colors[id]=(float3)(0.0f, 0.0f, 0.0f);
    }

    //if(id>180000 && id<180004){
    //    printf("colors[%06d]: [%06.2f %06.2f %06.2f]\n\r", id, colors[id].x, colors[id].y, colors[id].z);
    //    printf("hits[%06d]=\tt=%06.2f \tP=[%06.2f %06.2f %06.2f] \tN=[%06.2f %06.2f %06.2f]\n\r", id, hit.t, hit.P.x, hit.P.y, hit.P.z, hit.N.x, hit.N.y, hit.N.z);
    //}
    //printf("colors[%06d]: [%06.2f %06.2f %06.2f]\t", id, colors[id].x, colors[id].y, colors[id].z);
    //printf("hits[%06d]=\tt=%06.2f \tP=[%06.2f %06.2f %06.2f] \tN=[%06.2f %06.2f %06.2f]\n\r", id, hit.t, hit.P.x, hit.P.y, hit.P.z, hit.N.x, hit.N.y, hit.N.z);

    //printf("old P=[%06.2f %06.2f %06.2f] \tD=[%06.2f %06.2f %06.2f]\n\r", old_ray.P.x, old_ray.P.y, old_ray.P.z, old_ray.D.x, old_ray.D.y, old_ray.D.z);
    //printf("new P=[%06.2f %06.2f %06.2f] \tD=[%06.2f %06.2f %06.2f]\n\r", new_ray.P.x, new_ray.P.y, new_ray.P.z, new_ray.D.x, new_ray.D.y, new_ray.D.z);

    //printf("id=%02d \t%f %d\n\r", id, hit.mat.n, hit.mat.type);
    //printf("rays[%03d]=\tP=[%06.2f %06.2f %06.2f] \tD=[%06.2f %06.2f %06.2f]\n\r", id, rays[id].P.x, rays[id].P.y, rays[id].P.z, rays[id].D.x, rays[id].D.y, rays[id].D.z);
}

void kernel gen_ray(global Ray* rays, const Camera camera){
    int id=get_global_id(0);
    rays[id]=camera_get_ray(id, camera);
}