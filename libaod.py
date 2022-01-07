from pylab import *
import scipy.special as sp
from magic.libmagic import rderavg, thetaderavg

def gen_arr(lmax, l1,m1,mode='g'):

    idx = zeros([lmax+1,lmax+1])
    lArr = []
    mArr = []

    count = 0

    glm = []
    hlm = []

    for l in range(lmax+1):
        for m in range(l+1):

            if l in l1 and m in m1:
                if mode == 'g' or mode == 'gh':
                    glm.append(1.)
                else:
                    glm.append(0.)
                if mode == 'h' or mode == 'gh':
                    hlm.append(1.)
                else:
                    hlm.append(0.)
            else:
                glm.append(0.)
                hlm.append(0.)

            idx[l,m] = count
            lArr.append(l)
            mArr.append(m)
            
            count += 1

    glm  = array(glm)
    hlm  = array(hlm)
    lArr = array(lArr)
    mArr = array(mArr)
    idx = int32(idx)

    return glm, hlm, lArr, mArr, idx

def get_grid(r,theta,phi,nr,nt,np):
    r3D  = zeros([np,nt,nr])
    th3D = zeros([np,nt,nr])
    p3D  = zeros([np,nt,nr])
    for i in range(nr):
        r3D[...,i] = r[i]
    for j in range(nt):
        th3D[:,j,:] = theta[j]
    for k in range(np):
        p3D[k,...] = phi[k]
    s = r3D * sin(th3D)
    x = s * cos(p3D)
    y = s * sin(p3D)
    z = r3D * cos(th3D)
    return r3D,th3D,p3D, x,y,z, s

def getB_SPH(lmax,glm,hlm,idx,r3D,p3D,th3D,ro):
    
    Br = zeros_like(p3D)

    for l in range(lmax+1):
        for m in range(l+1):
            ylm = (-1)**m * sp.sph_harm(m, l, p3D, th3D)
            fac = (l+1)* (ro/r3D)**(l+2) * sqrt((4*pi)/(2*l+1))
            
            G = glm[idx[l,m]] * real(ylm)
            H = hlm[idx[l,m]] * imag(ylm) 

            Br +=   real(fac * (G + H))

    return Br

def getB_SPH_mer(r3D,p3D,th3D,ro):
    lmax = 2
    l1 = [1, 2]
    m1 = [0, 0]

    glm, hlm, lArr, mArr, idx = gen_arr(lmax, l1, m1)
    glm[idx[2,0]] = -74.6e-9
    glm[idx[1,0]] = -190e-9

    B = getB_SPH(lmax,glm,hlm,idx,r3D,p3D,th3D,ro)

    return B

def BCart2Sph(Bx,By,Bz,th3D,p3D):
    Br = (Bx * sin(th3D) * cos(p3D)) + (By * sin(th3D) * sin(p3D)) + (Bz * cos(th3D))
    Bt = (Bx * cos(th3D) * cos(p3D)) + (By * cos(th3D) * sin(p3D)) + (Bz * sin(th3D))
    Bp = -(Bx * sin(p3D)) + (By * cos(p3D))
    return Br, Bt, Bp

def curlVec(vr,vt,vp,r3D,th3D,eta,theta):

    curl_r = zeros_like(r3D)

    # Avoid the poles, d/dphi = 0

    curl_r[:,1:-1,:] = 1./(r3D[:,1:-1,:] * sin(th3D[:,1:-1,:])) \
                          * gradient(vp[:,1:-1,:]*sin(th3D[:,1:-1,:]), theta[1:-1], axis=1, edge_order=2)

    #curl_r[:,1:-1,:] = 1./(r3D[:,1:-1,:] * sin(th3D[:,1:-1,:])) * ( thetaderavg(vp[:,1:-1,:] * sin(th3D[:,1:-1,:])) )
    curl_t = 1./r3D * ( - rderavg(r3D * vp, eta = eta, spectral=True) )
    curl_p = 1./r3D * ( rderavg(r3D * vt, eta=eta, spectral=True) - thetaderavg(vr)  )

    return curl_r,curl_t,curl_p 

def BCart(x,y,z,m_dip, m_quad, th3D, p3D):

    r_mag = sqrt(x**2 + y**2 + z**2)

    fac_dip = m_dip/r_mag**5
    fac_quad = 1.5*m_quad/r_mag**7

    Bx = fac_dip * (3 * z * x)            + fac_quad * (5 * z**2 - r_mag**2) * x
    By = fac_dip * (3 * z * y)            + fac_quad * (5 * z**2 - r_mag**2) * y
    Bz = fac_dip * (3 * z * z - r_mag**2) + fac_quad * (5 * z**3 - 3 * r_mag**2 * z)

    Br, Bt, Bp = BCart2Sph(Bx,By,Bz,th3D,p3D)

    return Br

def getValf(x,y,z,m_dip, m_quad, th3D, p3D):
    
    r_mag = sqrt(x**2 + y**2 + z**2)

    fac_dip = m_dip/r_mag**5
    fac_quad = 1.5*m_quad/r_mag**7

    Bx = fac_dip * (3 * z * x)            + fac_quad * (5 * z**2 - r_mag**2) * x
    By = fac_dip * (3 * z * y)            + fac_quad * (5 * z**2 - r_mag**2) * y
    Bz = fac_dip * (3 * z * z - r_mag**2) + fac_quad * (5 * z**3 - 3 * r_mag**2 * z)

    B = sqrt(Bx**2 + By**2 + Bz**2)

    mu0 = 4*pi*1e-7

    rho = 10 * 1e6 * 1.67e-27

    Valf = B/sqrt(mu0*rho)

    return Valf

def VecPotential(x,y,z,m_dip,m_quad,th3D,p3D):

    r_mag = sqrt(x**2 + y**2 + z**2)

    fac_dip = m_dip/r_mag**3
    fac_quad = 1.5 * m_quad/r_mag**5

    Ax = fac_dip * (-y) + fac_quad * (- z* y)
    Ay = fac_dip *   x  + fac_quad * (z * x)
    Az = zeros_like(x)

    return Ax, Ay, Az

def BVecPot(x,y,z,m_dip, m_quad, r3D, th3D, p3D,eta,theta):

    Ax, Ay, Az = VecPotential(x,y,z,m_dip,m_quad,th3D,p3D)

    Ar,Atheta,Aphi = BCart2Sph(Ax,Ay,Az,th3D,p3D)   

    Br, Btheta, Bphi = curlVec(Ar,Atheta,Aphi,r3D,th3D,eta,theta)

    return Br

def fixClim(dat):

    datMax = max(abs(dat.max()),abs(dat.min()))
    clim(-datMax,datMax)

def merContour(r,theta,dat,levels=30,cmap='RdBu_r',colbar=True):
    theta2D, r2D = np.meshgrid(theta,r,indexing='ij')
    xx = r2D * np.sin(theta2D)
    yy = r2D * np.cos(theta2D)
    cont = contourf(xx,yy,dat,levels,cmap=cmap)
    axhline(y=0, ls=':', color = 'k')
    fixClim(dat)
    if colbar:
        cbar = colorbar(cont)
    
    for c in cont.collections:
        c.set_edgecolor("face")
    
    axis('equal')
    axis('off')
