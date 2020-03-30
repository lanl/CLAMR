inline real_t U_fullstep_version(
        real_t    deltaT,
        real_t    dr,
        real_t    U,
        real_t    F_plus,
        real_t    F_minus,
        real_t    G_plus,
        real_t    G_minus,
        real_t    wplusx_H,
        real_t    wminusx_H,
        real_t    wplusy_H,
        real_t    wminusy_H
) {

#if ( UPDATE_EQUATION_VERSION == 0 )
return U - (deltaT / dr)*(F_plus - F_minus + G_plus - G_minus) + - wminusx_H + wplusx_H  + - wminusy_H + wplusy_H;
#endif

#if ( UPDATE_EQUATION_VERSION == 1 )
return ((((U + (deltaT/dr)*F_minus) + (deltaT/dr)*G_minus) + wplusx_H) + wplusy_H) + (((-(deltaT/dr)*F_plus + -(deltaT/dr)*G_plus) + -wminusx_H) + -wminusy_H ); 
#endif


//more here


//If UPDATE_EQUATION_VERSION was not equal to any listed options, return default equation version
//Perhaps this should throw a compile warning or error?
return U - (deltaT / dr)*(F_plus - F_minus + G_plus - G_minus) + - wminusx_H + wplusx_H  + - wminusy_H + wplusy_H;

}
