""" Short MADX Codes """


def get_ptc_twiss_rdt(path, id):
    """ Create PTC universe, do ptc_twiss and output twissrdts, twiss, summary and nonlin tables """
    return (
        "ptc_create_universe;\n"
        "    ptc_create_layout, model=1, method=6, nst=5, exact=true, closed_layout=true;\n"
        "    ptc_setswitch, debuglevel=1, exact_mis=true, time=true, totalpath=false;\n"
        "    ptc_twiss, table=twiss, icase=5, no=4, closed_orbit=false, rmatrix, normal, writetmap, trackrdts;\n"
        "    write, table=twissrdt, file='%(OUTPUT_DIR)s/ptc_twiss.rdts.%(ID)s.dat';\n"
        "    write, table=twiss, file='%(OUTPUT_DIR)s/ptc_twiss.twiss.%(ID)s.dat';\n"
        "    write, table=ptc_twiss_summary, file='%(OUTPUT_DIR)s/ptc_twiss.summary.%(ID)s.dat';\n"
        "    write, table=nonlin, file='%(OUTPUT_DIR)s/ptc_twiss.nonlin.%(ID)s.dat';\n"
        "ptc_end;\n") % {"OUTPUT_DIR": path, "ID": id}


def get_ptc_amplitude_detuning(path, id):
    """ Create PTC universe, do ptc_normal for amplitude detuning """
    return (
        "ptc_create_universe;\n"
        "    ptc_create_layout,model=3,method=6,nst=3,resplit,thin=0.0005,xbend=0.0005;\n"
        "    ptc_align;\n"
        "!    ptc_setswitch, fringe=True;\n"
        "    select_ptc_normal,  q1=0, q2=0;\n"
        "    select_ptc_normal, dq1=1,dq2=1;\n"
        "    select_ptc_normal, dq1=2,dq2=2;\n"
        "    select_ptc_normal, dq1=3,dq2=3;\n"
        "    select_ptc_normal, anhx=1,0,0; ! dQx/dex\n"
        "    select_ptc_normal, anhy=0,1,0; ! dQy/dey\n"
        "    select_ptc_normal, anhx=0,1,0;\n"
        "    select_ptc_normal, anhy=1,0,0;\n"
        "    select_ptc_normal, anhx=2,0,0; ! d2Qx/dex^2\n"
        "    select_ptc_normal, anhx=1,1,0;\n"
        "    select_ptc_normal, anhx=0,2,0;\n"
        "    select_ptc_normal, anhy=0,2,0; ! d2Qy/dey^2\n"
        "    select_ptc_normal, anhy=1,1,0; ! d2Qy/deydex\n"
        "    select_ptc_normal, anhy=2,0,0;\n"
        "    ptc_normal,closed_orbit,normal,icase=5,no=5;\n"
        "    write, table=normal_results,file='%(OUTPUT_DIR)s/ptc_normal.ampdet.%(ID)s.dat';\n"
        "ptc_end;\n") % {"OUTPUT_DIR": path, "ID": id}


# Script Mode ##################################################################


if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))
