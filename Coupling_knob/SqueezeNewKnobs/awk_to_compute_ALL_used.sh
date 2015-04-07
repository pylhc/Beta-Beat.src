

awk '/dq/{print FILENAME, a/$3**2*1e-6;a=0}/kq/{a=a+$3**2}' used_c-imag_b2.IP?.dat > ti

awk '/dq/{print FILENAME, a/$3**2*1e-6;a=0}/kq/{a=a+$3**2}' used_c-real_b2.IP?.dat > tr


paste tr ti > tt

awk '{print $2+$4, $0}' tt > ALL_used_Strengths.b2.dat


rm tr ti tt
