
typedef struct{
    float3 kd,ks,emission,F0;
    float n, shininess, glossiness;
    int type;
} Material;

typedef struct {
    float3 P,D;
} Ray;

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

typedef struct {
    float3 r1,r2,r3,N;
    Material mat;
} Triangle;
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
        //printf("hits[%03d]=\tt=%06.2f \tP=[%06.2f %06.2f %06.2f] \tN=[%06.2f %06.2f %06.2f]\n", i, hit.t, hit.P.x, hit.P.y, hit.P.z, hit.N.x, hit.N.y, hit.N.z);
        if(hit.t>0 && (best_hit.t<0 || hit.t<best_hit.t))
            best_hit=hit;
    }
    return best_hit;
}

void kernel trace_ray(global const Triangle* tris, const int tris_size, global const Ray* rays, global Hit* hits){
    int id=get_global_id(0);
    Hit hit=first_intersect(tris, tris_size, rays[id]);
    hits[id]=hit;
}