
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
                return cons_Hit(t,p,N);
            }
        }
    }
    return init_Hit();
}

void kernel trace_ray(global const Triangle* tris, const int tris_size, global const Ray* rays, global Hit* hits){
    //Triangle tri=tris[0];
    //printf("tris[%d]= r1=[%2.2f %2.2f %2.2f] r2=[%2.2f %2.2f %2.2f] r3=[%2.2f %2.2f %2.2f] N=[%2.2f %2.2f %2.2f]\t", tris_size, tri.r1.x, tri.r1.y, tri.r1.z, tri.r2.x, tri.r2.y, tri.r2.z, tri.r3.x, tri.r3.y, tri.r3.z, tri.N.x, tri.N.y, tri.N.z);
    //printf("tris_size= %d\t",tris_size);

    //Ray ray=rays[id];
    //printf("rays[%d]= P=[%2.2f %2.2f %2.2f] D=[%2.2f %2.2f %2.2f]\n", id, ray.P.x, ray.P.y, ray.P.z, ray.D.x, ray.D.y, ray.D.z);
    
    //Hit hit=hits[get_global_id(0)];
    //printf("hits[0]= t=%2.2f P=[%2.2f %2.2f %2.2f] N=[%2.2f %2.2f %2.2f]\n", hit.t, hit.P.x, hit.P.y, hit.P.z, hit.N.x, hit.N.y, hit.N.z);


    int id=get_global_id(0);
    Hit bestHit=init_Hit();
    for(int i=0; i<tris_size; ++i){
        Hit hit=triangle_intersect(tris[i], rays[id]);
        //printf("hits[%03d]=\tt=%06.2f \tP=[%06.2f %06.2f %06.2f] \tN=[%06.2f %06.2f %06.2f]\n", i, hit.t, hit.P.x, hit.P.y, hit.P.z, hit.N.x, hit.N.y, hit.N.z);
        if(hit.t>0 && (bestHit.t<0 || hit.t<bestHit.t))
            bestHit=hit;
    }
    hits[id]=bestHit;
}