#include <bits/stdc++.h>
using namespace std;
using u64=uint64_t;using u128=__uint128_t;

struct Field {
    int p; vector<int> inv;
    Field(int p_=2,int n=0):p(p_){ if(n) build_inv(n); }
    int norm(long long x)const{x%=p;if(x<0)x+=p;return (int)x;}
    int add(int a,int b)const{int s=a+b;if(s>=p)s-=p;return s;}
    int sub(int a,int b)const{int s=a-b;if(s<0)s+=p;return s;}
    int mul(long long a,long long b)const{return int(a*b%p);}
    int powmod(int a,long long e)const{long long r=1,b=a;while(e){if(e&1)r=r*b%p;b=b*b%p;e>>=1;}return (int)r;}
    int inverse(int a)const{return powmod(a,p-2);}
    void build_inv(int n){inv.assign(n+1,0);for(int i=1;i<=n;i++){int a=i%p;if(a)inv[i]=inverse(a);}}
};

template<int MOD,int ROOT>
struct NTT {
    static int add(int a,int b){long long s=(long long)a+b;if(s>=MOD)s-=MOD;return (int)s;}
    static int sub(int a,int b){long long s=(long long)a-b;if(s<0)s+=MOD;return (int)s;}
    static int mul(long long a,long long b){return int(a*b%MOD);}
    static int pw(int a,long long e){long long r=1,b=a;while(e){if(e&1)r=r*b%MOD;b=b*b%MOD;e>>=1;}return (int)r;}
    static void transform(vector<int>&a,bool invert){
        int n=a.size();
        for(int i=1,j=0;i<n;i++){int bit=n>>1;for(;j&bit;bit>>=1)j^=bit;j^=bit;if(i<j)swap(a[i],a[j]);}
        for(int len=2;len<=n;len<<=1){int wl=pw(ROOT,(MOD-1)/len);if(invert)wl=pw(wl,MOD-2);for(int i=0;i<n;i+=len){long long w=1;int h=len>>1;for(int j=0;j<h;j++){int u=a[i+j],v=mul(a[i+j+h],w);a[i+j]=add(u,v);a[i+j+h]=sub(u,v);w=w*wl%MOD;}}}
        if(invert){int ni=pw(n,MOD-2);for(int&x:a)x=mul(x,ni);}
    }
};
using N1=NTT<2013265921,31>; // 15*2^27+1
using N2=NTT<469762049,3>;   // 7*2^26+1
static constexpr u64 M1=2013265921ULL,M2=469762049ULL;

u64 crt2(int a,int b){
    static const u64 invM1modM2=[](){
        long long aa=M1%M2, e=M2-2, r=1, x=aa;while(e){if(e&1)r=r*x%M2;x=x*x%M2;e>>=1;}return (u64)r;
    }();
    long long diff=(long long)b-(long long)(a%M2);diff%= (long long)M2;if(diff<0)diff+=M2;
    u64 t=(u128)(u64)diff*invM1modM2%M2;
    return (u64)a+M1*t;
}

vector<int> mobius_table(int N){vector<int>mu(N+1),lp(N+1),ps;mu[1]=1;for(int i=2;i<=N;i++){if(!lp[i]){lp[i]=i;ps.push_back(i);mu[i]=-1;}for(int p:ps){if(p>lp[i]||1LL*i*p>N)break;lp[i*p]=p;mu[i*p]=(i%p? -mu[i]:0);}}return mu;}
vector<vector<int>> divisor_table(int N){vector<vector<int>>d(N+1);for(int x=1;x<=N;x++)for(int y=x;y<=N;y+=x)d[y].push_back(x);return d;}

struct LowTable {
    int G,W,S,p;
    Field F;
    vector<int> inv,mu,Bern,fact,ifact;
    vector<vector<int>> divs;
    vector<vector<int>> R; // [w][tri]
    vector<vector<int>> H; // [h][w]
    vector<vector<int>> E; // [w][(g-1)*(S+1)+s], output coefficient by w degree

    LowTable(int G_,int W_,int S_,int p_):G(G_),W(W_),S(S_),p(p_),F(p_,max(2*G_,S_)){
        if(p<=max(G+2,S)) throw runtime_error("p must exceed G+2 and S");
        inv=F.inv;mu=mobius_table(2*G);divs=divisor_table(2*G);
        fact.resize(G+2);ifact.resize(fact.size());fact[0]=1;for(int i=1;i<(int)fact.size();i++)fact[i]=F.mul(fact[i-1],i);ifact.back()=F.inverse(fact.back());for(int i=(int)fact.size()-1;i>=1;i--)ifact[i-1]=F.mul(ifact[i],i);
        build_bernoulli(G);
    }
    int C(int n,int k)const{if(k<0||k>n)return 0;return F.mul(fact[n],F.mul(ifact[k],ifact[n-k]));}
    int tri(int m,int j)const{return m*(m+1)/2+j;}
    void build_bernoulli(int N){Bern.assign(N+1,0);Bern[0]=1;for(int m=1;m<=N;m++){int sm=0;for(int k=0;k<m;k++)sm=F.add(sm,F.mul(C(m+1,k),Bern[k]));Bern[m]=F.sub(0,F.mul(sm,inv[m+1]));}}

    void build_R(){
        int base=2*G+1,Emax=G*base+G,N=1;while(N<=2*Emax)N<<=1;if(N>(1<<26))throw runtime_error("R NTT too large for selected NTT primes");
        cerr<<"R N="<<N<<"\n";
        vector<vector<int>> K(W+1,vector<int>(N));
        auto enc=[&](int m,int j){return m*base+j;};
        for(int m=1;m<=G;m++){K[0][enc(m,0)]=F.sub(K[0][enc(m,0)],inv[m]);if(W>=1)K[1][enc(m,0)]=F.add(K[1][enc(m,0)],inv[m]);}
        vector<int> Rs(2*G+1),logRs(G+1),fpow(2*G+1);
        for(int ell=2;ell<=2*G;ell++){if(ell%p==0)continue;
            vector<int> wp(W+1);int iell=inv[ell];for(int d:divs[ell])if(d<=W&&mu[ell/d])wp[d]=F.add(wp[d],F.mul(F.norm(-mu[ell/d]),iell));bool nz=false;for(int b=1;b<=W;b++)nz|=wp[b]!=0;if(!nz)continue;
            int c=F.mul(F.norm(mu[ell]),iell);
            fill(Rs.begin(),Rs.end(),0);Rs[0]=1;vector<pair<int,int>>supp;for(int d:divs[ell])if(d<ell&&mu[ell/d]){int e=ell-d;if(e<=2*G)Rs[e]=F.add(Rs[e],F.norm(mu[ell/d]));}for(int e=1;e<=2*G;e++)if(Rs[e])supp.push_back({e,Rs[e]});
            fill(logRs.begin(),logRs.end(),0);for(int n=1;n<=G;n++){int sm=0;for(auto[a,av]:supp){if(a>=n)break;int k=n-a;sm=F.add(sm,F.mul(F.mul(k,logRs[k]),av));}int ln=F.sub(Rs[n],F.mul(sm,inv[n]));logRs[n]=ln;int h=ln;if(n%ell==0)h=F.sub(h,inv[n/ell]);if(h)for(int b=1;b<=W;b++)if(wp[b])K[b][enc(n,0)]=F.add(K[b][enc(n,0)],F.mul(wp[b],h));}
            vector<vector<int>>wpow(W+1,vector<int>(W+1));wpow[0][0]=1;for(int a=1;a<=W;a++)for(int i=0;i<=W;i++)if(wpow[a-1][i])for(int j=1;i+j<=W;j++)if(wp[j])wpow[a][i+j]=F.add(wpow[a][i+j],F.mul(wpow[a-1][i],wp[j]));
            vector<int>cpow(G+W+2);cpow[0]=1;for(int j=1;j<(int)cpow.size();j++)cpow[j]=F.mul(cpow[j-1],c);
            int ellpow=1,smax=2*G/ell;
            for(int ss=1;ss<=smax;ss++){
                ellpow=F.mul(ellpow,ell);int nf=2*G-ell*ss;fill(fpow.begin(),fpow.begin()+nf+1,0);fpow[0]=1;
                for(int n=1;n<=nf;n++){int sm=0;for(auto[a,av]:supp){if(a>n)break;int wt=F.norm(n+1LL*(ss-1)*a);sm=F.add(sm,F.mul(F.mul(wt,av),fpow[n-a]));}fpow[n]=F.sub(0,F.mul(sm,inv[n]));}
                auto proc=[&](int k,int coef){if(k<1||coef==0)return;int amin=max(1,k-G),amax=min(k,W);for(int a=amin;a<=amax;a++){int j=k-a;int sc=F.mul(coef,F.mul(C(k,a),cpow[j]));if(!sc)continue;int lo=max(0,2*j-ell*ss),hi=min(nf,G+j-ell*ss);if(lo>hi)continue;for(int b=a;b<=W;b++)if(wpow[a][b]){int sb=F.mul(sc,wpow[a][b]);for(int nn=lo;nn<=hi;nn++)if(fpow[nn]){int m=ell*ss+nn-j;int v=F.mul(sb,F.mul(ellpow,fpow[nn]));K[b][enc(m,j)]=F.add(K[b][enc(m,j)],v);}}}};
                proc(ss,inv[2*ss]);proc(ss+1,F.sub(0,F.mul(inv[ss],inv[ss+1])));for(int rb=2;rb<=ss;rb+=2){int k=ss-rb+1;int co=F.sub(0,F.mul(Bern[rb],F.mul(C(ss-1,k),F.mul(inv[rb],inv[rb-1]))));proc(k,co);}
            }
        }
        // arbitrary-modulus w-exponential via two exact integer NTT convolutions
        u128 bound=(u128)(W*(W+1)/2)*(u128)((G+1)*(G+2)/2)*(u128)(p-1)*(p-1);
        if(bound >= (u128)M1*M2) throw runtime_error("CRT bound insufficient in R");
        vector<vector<int>>K1(W+1),K2(W+1),R1(W+1),R2(W+1);for(int b=1;b<=W;b++){K1[b]=K[b];K2[b]=K[b];N1::transform(K1[b],false);N2::transform(K2[b],false);}R.assign(W+1,vector<int>((G+1)*(G+2)/2));
        R1[0].assign(N,0);R2[0].assign(N,0);R1[0][0]=R2[0][0]=1;R1[0][enc(1,0)]=R2[0][enc(1,0)]=p-1;R[0][tri(0,0)]=1;R[0][tri(1,0)]=p-1;N1::transform(R1[0],false);N2::transform(R2[0],false);
        vector<int>a1(N),a2(N),coef1(N),coef2(N);
        for(int a=1;a<=W;a++){
            fill(a1.begin(),a1.end(),0);fill(a2.begin(),a2.end(),0);for(int b=1;b<=a;b++){for(int i=0;i<N;i++){a1[i]=N1::add(a1[i],N1::mul(b,N1::mul(K1[b][i],R1[a-b][i])));a2[i]=N2::add(a2[i],N2::mul(b,N2::mul(K2[b][i],R2[a-b][i])));}}N1::transform(a1,true);N2::transform(a2,true);fill(coef1.begin(),coef1.end(),0);fill(coef2.begin(),coef2.end(),0);for(int m=0;m<=G;m++)for(int j=0;j<=m;j++){int idx=enc(m,j);int v=(int)(crt2(a1[idx],a2[idx])%p);v=F.mul(v,inv[a]);R[a][tri(m,j)]=v;coef1[idx]=coef2[idx]=v;}if(a<W){R1[a]=coef1;R2[a]=coef2;N1::transform(R1[a],false);N2::transform(R2[a],false);}}
    }

    void build_H(){
        H.assign(G,vector<int>(W+1));H[0][0]=1;vector<vector<int>>d(G,vector<int>(W+1));
        for(int q=1;q<G;q++)for(int a=1;a<=W&&a<=q+1;a++){int idx=q+1-a;int bx=(idx==1?F.mul(1,inv[2]):Bern[idx]);int val=F.mul((a&1)?F.norm(-1):1,F.mul(C(q+1,a),bx));val=F.sub(0,F.mul(val,F.mul(inv[q],inv[q+1])));d[q][a]=val;}
        for(int q=1;q<G;q++){vector<int>acc(W+1);for(int rr=1;rr<=q;rr++)for(int a=1;a<=W;a++)if(d[rr][a])for(int b=0;a+b<=W;b++)if(H[q-rr][b])acc[a+b]=F.add(acc[a+b],F.mul(rr,F.mul(d[rr][a],H[q-rr][b])));for(int a=0;a<=W;a++)H[q][a]=F.mul(acc[a],inv[q]);}
    }

    vector<vector<int>> build_A_coeff(){
        int width=S+1;vector<vector<int>>A(W+1,vector<int>(G*width));vector<int>poly(W+1),npoly(W+1),prod(W+1);
        for(int h=0;h<G;h++){
            fill(poly.begin(),poly.end(),0);poly[0]=1;
            for(int n=0;n<=S;n++){
                fill(prod.begin(),prod.end(),0);for(int a=0;a<=W;a++)if(H[h][a])for(int b=0;a+b<=W;b++)if(poly[b])prod[a+b]=F.add(prod[a+b],F.mul(H[h][a],poly[b]));for(int a=0;a<=W;a++)A[a][h*width+n]=prod[a];
                if(n<S){fill(npoly.begin(),npoly.end(),0);int c=F.norm(h-1+n);for(int a=0;a<=W;a++)if(poly[a]){npoly[a]=F.add(npoly[a],F.mul(c,poly[a]));if(a+1<=W)npoly[a+1]=F.add(npoly[a+1],poly[a]);}int iv=inv[n+1];for(int a=0;a<=W;a++)npoly[a]=F.mul(npoly[a],iv);poly.swap(npoly);}
            }
        }
        return A;
    }

    void multiply_AR(){
        auto A=build_A_coeff();int base=S+G+1;int maxR=(G-1)*base+G,maxA=(G-1)*base+S,N=1;while(N<=maxR+maxA)N<<=1;if(N>(1<<26))throw runtime_error("A*R NTT too large for selected NTT primes");cerr<<"AR N="<<N<<"\n";
        u128 ar_bound=(u128)(W+1)*(u128)((G+1)*(G+2)/2)*(u128)(p-1)*(p-1);
        if(ar_bound >= (u128)M1*M2) throw runtime_error("CRT bound insufficient in A*R");
        size_t outsz=(size_t)G*(S+1);vector<vector<int>>res1(W+1,vector<int>(outsz));E.assign(W+1,vector<int>(outsz));
        auto process1=[&](){vector<vector<int>>Ah(W+1,vector<int>(N)),Rh(W+1,vector<int>(N));for(int b=0;b<=W;b++){for(int h=0;h<G;h++)for(int n=0;n<=S;n++)Ah[b][h*base+n]=A[b][h*(S+1)+n];for(int m=0;m<G;m++)for(int j=0;j<=m;j++)Rh[b][m*base+j]=R[b][tri(m,j)];N1::transform(Ah[b],false);N1::transform(Rh[b],false);}vector<int>out(N);for(int a=0;a<=W;a++){fill(out.begin(),out.end(),0);for(int b=0;b<=a;b++)for(int i=0;i<N;i++)out[i]=N1::add(out[i],N1::mul(Ah[b][i],Rh[a-b][i]));N1::transform(out,true);for(int g=1;g<=G;g++)for(int s=0;s<=S;s++)res1[a][(g-1)*(S+1)+s]=out[(g-1)*base+s];}};
        process1();
        // Free A? retained, but memory acceptable. Process second modulus.
        vector<vector<int>>Ah(W+1,vector<int>(N)),Rh(W+1,vector<int>(N));for(int b=0;b<=W;b++){for(int h=0;h<G;h++)for(int n=0;n<=S;n++)Ah[b][h*base+n]=A[b][h*(S+1)+n];for(int m=0;m<G;m++)for(int j=0;j<=m;j++)Rh[b][m*base+j]=R[b][tri(m,j)];N2::transform(Ah[b],false);N2::transform(Rh[b],false);}vector<int>out(N);for(int a=0;a<=W;a++){fill(out.begin(),out.end(),0);for(int b=0;b<=a;b++)for(int i=0;i<N;i++)out[i]=N2::add(out[i],N2::mul(Ah[b][i],Rh[a-b][i]));N2::transform(out,true);for(int g=1;g<=G;g++)for(int s=0;s<=S;s++){size_t pos=(g-1)*(S+1)+s;int v=(int)(crt2(res1[a][pos],out[(g-1)*base+s])%p);E[a][pos]=v;}}
        // subtract the constant 1 from e^S-1
        E[0][0]=F.sub(E[0][0],1);
    }

    void run(){auto t=chrono::steady_clock::now();build_R();cerr<<"R sec "<<chrono::duration<double>(chrono::steady_clock::now()-t).count()<<"\n";t=chrono::steady_clock::now();build_H();cerr<<"H sec "<<chrono::duration<double>(chrono::steady_clock::now()-t).count()<<"\n";t=chrono::steady_clock::now();multiply_AR();cerr<<"AR sec "<<chrono::duration<double>(chrono::steady_clock::now()-t).count()<<"\n";}
    int coeff(int r,int g,int s)const{int v=0;size_t pos=(g-1)*(S+1)+s;for(int a=0;a<=r;a++)v=F.add(v,E[a][pos]);return v;}
};

int main(int argc,char**argv){if(argc<5){cerr<<"usage: modular_low_truncation_table_large G W S p [output.bin]\n";return 2;}int G=stoi(argv[1]),W=stoi(argv[2]),S=stoi(argv[3]),p=stoi(argv[4]);auto t=chrono::steady_clock::now();LowTable T(G,W,S,p);T.run();cerr<<"total "<<chrono::duration<double>(chrono::steady_clock::now()-t).count()<<"\n";u64 chk=0;for(int r=0;r<=W;r++)for(int g=1;g<=G;g++)for(int s=0;s<=S;s++)chk=(chk*1000003+T.coeff(r,g,s))%1000000007;cout<<"PASS checksum="<<chk<<"\n";if(argc>=6){ofstream f(argv[5],ios::binary);uint32_t hdr[4]={(uint32_t)G,(uint32_t)W,(uint32_t)S,(uint32_t)p};f.write((char*)hdr,sizeof(hdr));for(int r=0;r<=W;r++)for(int g=1;g<=G;g++)for(int s=0;s<=S;s++){uint16_t v=T.coeff(r,g,s);f.write((char*)&v,sizeof(v));}}}
