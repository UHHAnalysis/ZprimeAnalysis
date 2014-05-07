
### Filter definitions ###

def narrow_resonances(hname):
    # Accept anything that there is neither of the signals
    if 'rsg' not in hname and 'zp' not in hname :
        return True
    # Reject RS gluons as signal
    elif 'rsg' in hname:
        return False
    # Process signal name
    pname = hname.split('__')[1]
    # reject wide reonances
    if 'w10p' in pname:
        return False
    # Accept only a few mass points (no interpolation)
    mass = pname.split('w')[0].split('zp')[1]
    #mass_whitelist = ['500', '750','1000','1250','1500','2000','3000','4000']
    return float(mass) <= 3000


def wide_resonances(hname):
    # Accept anything that there is neither of the signals
    if 'rsg' not in hname and 'zp' not in hname :
        return True
    # Reject RS gluons as signal
    elif 'rsg' in hname:
        return False
    # Process signal name
    pname = hname.split('__')[1]
    # reject wide reonances
    if 'w1p' in pname:
        return False
    # Accept only a few mass points (no interpolation)
    mass = pname.split('w')[0].split('zp')[1]
    #mass_whitelist = ['500', '750','1000','1250','1500','2000','3000','4000']
    return float(mass) <= 3000


def rsg_resonances(hname):
    # Accept anything that there is neither of the signals
    if 'rsg' not in hname and 'zp' not in hname :
        return True
    # Reject zp as signal
    elif 'zp' in hname:
        return False
    # Process signal name
    pname = hname.split('__')[1]
    # Accept only a few mass points (no interpolation)
    mass = pname[3:]
    #mass_whitelist = ['1000','1500','2000','2500','3000','3500','4000']
    return float(mass) <= 3000


def build_boosted_semileptonic_model(files, filter, signal, mcstat, eflag=False, muflag=False):
    """ Semileptonic high mass model"""
    model = build_model_from_rootfile(files, filter, include_mc_uncertainties = mcstat)
    model.fill_histogram_zerobins()
    model.set_signal_processes(signal)
    for p in model.processes:
        model.add_lognormal_uncertainty('lumi', math.log(1.026), p)

    model.add_lognormal_uncertainty('zj_rate', math.log(2.0), 'zlight')
    model.add_lognormal_uncertainty('wj_rate', math.log(1.5), 'wlight')
    model.add_lognormal_uncertainty('wj_rate', math.log(1.5), 'wb')
    model.add_lognormal_uncertainty('wj_rate', math.log(1.5), 'wc')
    model.add_lognormal_uncertainty('wb_rate', math.log(1.87), 'wb')
    model.add_lognormal_uncertainty('wc_rate', math.log(1.87), 'wc')
    model.add_lognormal_uncertainty('ttbar_rate', math.log(1.15), 'ttbar')
    model.add_lognormal_uncertainty('st_rate', math.log(1.5), 'singletop')
    model.add_lognormal_uncertainty('diboson_rate', math.log(1.5), 'diboson')

    if muflag:
        for obs in ['mu_0top0btag_mttbar','mu_1top0btag_mttbar']:
            for proc in ('wc', 'wb'):
                model.add_asymmetric_lognormal_uncertainty('scale_vjets', -math.log(1.577), math.log(0.710), proc, obs)
                model.add_asymmetric_lognormal_uncertainty('matching_vjets', -math.log(1.104), math.log(1.052), proc, obs)
        for obs in ['mu_0top1btag_mttbar','mu_0top2btag_mttbar','mu_1top1btag_mttbar','mu_1top2btag_mttbar']:
            for proc in ('wc', 'wb', 'wlight'):
                model.add_asymmetric_lognormal_uncertainty('scale_vjets', -math.log(1.577), math.log(0.710), proc, obs)
                model.add_asymmetric_lognormal_uncertainty('matching_vjets', -math.log(1.104), math.log(1.052), proc, obs)

    if eflag:
        #For categories with low statistics, use flat uncertainties instead of shape
        #  Template for the following lines:
        #    model.add_asymmetric_lognormal_uncertainty(sys, -math.log( _PLUS_ ), math.log( _MINUS_ ), proc, obs)
        #  Where _PLUS_ = proc_sys_plus.Integral() / proc.Integral()  after the kinematic selection, before Chi2
        #  Where _MINUS_ = proc_sys_minus.Integral() / proc.Integral()  after the kinematic selection, before Chi2
        for obs in ['el_1top0btag_mttbar','el_1top1btag_mttbar','el_1top2btag_mttbar']:
            model.add_asymmetric_lognormal_uncertainty('bmistag', -math.log(0.99828), math.log(1.00173), 'diboson', obs)
            model.add_asymmetric_lognormal_uncertainty('btageff', -math.log(1.00020), math.log(0.99979), 'diboson', obs)
            model.add_asymmetric_lognormal_uncertainty('elesf'  , -math.log(1.00049), math.log(0.99951), 'diboson', obs)
            model.add_asymmetric_lognormal_uncertainty('jec'    , -math.log(1.03533), math.log(0.95318), 'diboson', obs)
            model.add_asymmetric_lognormal_uncertainty('jer'    , -math.log(1.01404), math.log(0.98601), 'diboson', obs)
            model.add_asymmetric_lognormal_uncertainty('pileup' , -math.log(1.00604), math.log(0.99374), 'diboson', obs)
            model.add_asymmetric_lognormal_uncertainty('topmistag', -math.log(1.00007), math.log(0.99993), 'diboson', obs)
            model.add_asymmetric_lognormal_uncertainty('toptageff', -math.log(1.00004), math.log(0.99996), 'diboson', obs)

            model.add_asymmetric_lognormal_uncertainty('bmistag', -math.log(1.00023), math.log(0.99977), 'singletop', obs)
            model.add_asymmetric_lognormal_uncertainty('btageff', -math.log(0.99951), math.log(1.00049), 'singletop', obs)
            model.add_asymmetric_lognormal_uncertainty('elesf'  , -math.log(1.00044), math.log(0.99956), 'singletop', obs)
            model.add_asymmetric_lognormal_uncertainty('jec'    , -math.log(1.02887), math.log(0.95917), 'singletop', obs)
            model.add_asymmetric_lognormal_uncertainty('jer'    , -math.log(1.01413), math.log(0.98887), 'singletop', obs)
            model.add_asymmetric_lognormal_uncertainty('pileup' , -math.log(1.00121), math.log(0.99904), 'singletop', obs)
            model.add_asymmetric_lognormal_uncertainty('topmistag', -math.log(0.99995), math.log(1.00005), 'singletop', obs)
            model.add_asymmetric_lognormal_uncertainty('toptageff', -math.log(1.00036), math.log(0.99964), 'singletop', obs)

            model.add_asymmetric_lognormal_uncertainty('bmistag', -math.log(1.00036), math.log(0.99963), 'wb', obs)
            model.add_asymmetric_lognormal_uncertainty('btageff', -math.log(0.99718), math.log(1.00277), 'wb', obs)
            model.add_asymmetric_lognormal_uncertainty('elesf'  , -math.log(1.00051), math.log(0.99949), 'wb', obs)
            model.add_asymmetric_lognormal_uncertainty('jec'    , -math.log(1.03401), math.log(0.95899), 'wb', obs)
            model.add_asymmetric_lognormal_uncertainty('jer'    , -math.log(1.01522), math.log(0.99024), 'wb', obs)
            model.add_asymmetric_lognormal_uncertainty('pileup' , -math.log(1.00922), math.log(0.98990), 'wb', obs)
            model.add_asymmetric_lognormal_uncertainty('topmistag', -math.log(1.00031), math.log(0.99969), 'wb', obs)
            #model.add_asymmetric_lognormal_uncertainty('toptageff', -math.log(1.0), math.log(1.0), 'wb', obs)
            model.add_asymmetric_lognormal_uncertainty('scale_vjets', -math.log(2.317105), math.log(0.50625), 'wb', obs)
            model.add_asymmetric_lognormal_uncertainty('matching_vjets', -math.log(0.90346), math.log(0.64922), 'wb', obs)

            model.add_asymmetric_lognormal_uncertainty('bmistag', -math.log(0.99854), math.log(1.00146), 'wc', obs)
            model.add_asymmetric_lognormal_uncertainty('btageff', -math.log(0.99908), math.log(1.00093), 'wc', obs)
            model.add_asymmetric_lognormal_uncertainty('elesf'  , -math.log(1.00048), math.log(0.99951), 'wc', obs)
            model.add_asymmetric_lognormal_uncertainty('jec'    , -math.log(1.03686), math.log(0.95038), 'wc', obs)
            model.add_asymmetric_lognormal_uncertainty('jer'    , -math.log(1.01454), math.log(0.98724), 'wc', obs)
            model.add_asymmetric_lognormal_uncertainty('pileup' , -math.log(1.00449), math.log(0.99578), 'wc', obs)
            model.add_asymmetric_lognormal_uncertainty('topmistag', -math.log(1.00007), math.log(0.99993), 'wc', obs)
            #model.add_asymmetric_lognormal_uncertainty('toptageff', -math.log(1.0), math.log(1.0), 'wb', obs)
            model.add_asymmetric_lognormal_uncertainty('scale_vjets', -math.log(2.35657), math.log(0.49810), 'wc', obs)
            model.add_asymmetric_lognormal_uncertainty('matching_vjets', -math.log(0.87614), math.log(0.67215), 'wc', obs)

            model.add_asymmetric_lognormal_uncertainty('bmistag', -math.log(0.99840), math.log(1.00162), 'wlight', obs)
            model.add_asymmetric_lognormal_uncertainty('btageff', -math.log(0.99992), math.log(1.00008), 'wlight', obs)
            model.add_asymmetric_lognormal_uncertainty('elesf'  , -math.log(1.00051), math.log(0.99949), 'wlight', obs)
            model.add_asymmetric_lognormal_uncertainty('jec'    , -math.log(1.03593), math.log(0.95568), 'wlight', obs)
            model.add_asymmetric_lognormal_uncertainty('jer'    , -math.log(1.01414), math.log(0.98731), 'wlight', obs)
            model.add_asymmetric_lognormal_uncertainty('pileup' , -math.log(1.00495), math.log(0.99499), 'wlight', obs)
            model.add_asymmetric_lognormal_uncertainty('topmistag', -math.log(0.99998), math.log(1.00002), 'wlight', obs)
            #model.add_asymmetric_lognormal_uncertainty('toptageff', -math.log(1.0), math.log(1.0), 'wlight', obs)
            model.add_asymmetric_lognormal_uncertainty('scale_vjets', -math.log(2.32356), math.log(0.50250), 'wlight', obs)
            model.add_asymmetric_lognormal_uncertainty('matching_vjets', -math.log(0.89011), math.log(0.69660), 'wlight', obs)

            model.add_asymmetric_lognormal_uncertainty('bmistag', -math.log(1.00019), math.log(0.99979), 'zlight', obs)
            model.add_asymmetric_lognormal_uncertainty('btageff', -math.log(0.99975), math.log(1.00025), 'zlight', obs)
            model.add_asymmetric_lognormal_uncertainty('elesf'  , -math.log(1.00050), math.log(0.99950), 'zlight', obs)
            model.add_asymmetric_lognormal_uncertainty('jec'    , -math.log(1.05828), math.log(0.94517), 'zlight', obs)
            model.add_asymmetric_lognormal_uncertainty('jer'    , -math.log(1.05491), math.log(0.96224), 'zlight', obs)
            model.add_asymmetric_lognormal_uncertainty('pileup' , -math.log(1.02015), math.log(0.98078), 'zlight', obs)
            model.add_asymmetric_lognormal_uncertainty('topmistag', -math.log(0.99995), math.log(1.00005), 'zlight', obs)
            #model.add_asymmetric_lognormal_uncertainty('toptageff', -math.log(1.0), math.log(1.0), 'zlight', obs)

        for obs in ['el_1top0btag_mttbar','el_1top1btag_mttbar','el_1top2btag_mttbar','el_0top0btag_mttbar','el_0top1btag_mttbar','el_0top2btag_mttbar']:
            model.add_asymmetric_lognormal_uncertainty('scale_vjets', -math.log(0.87969), math.log(1.21499), 'zlight', obs)
            model.add_asymmetric_lognormal_uncertainty('matching_vjets', -math.log(1.12298), math.log(1.67109), 'zlight', obs)
            for proc in model.processes:
                model.add_lognormal_uncertainty('eltrig_rate', math.log(1.01), p)

    return model


import exceptions


def build_model(type, jet1 = None, mcstat = True):

    model = None

    if type == 'narrow_resonances_muon':
        model = build_boosted_semileptonic_model(
            ['theta_input_muon_toptag_rebinned.root'],
            narrow_resonances,
            'zp*',
            mcstat,
            muflag = True
        )

    elif type == 'wide_resonances_muon':

        model = build_boosted_semileptonic_model(
            ['theta_input_muon_toptag_rebinned.root'],
            wide_resonances,
            'zp*',
            mcstat,
            muflag = True
        )

    elif type == 'rsg_resonances_muon':

        model = build_boosted_semileptonic_model(
            ['theta_input_muon_toptag_rebinned.root'],
            rsg_resonances,
            'rsg*',
            mcstat,
            muflag = True
        )

    elif type == 'narrow_resonances_electron':

        model = build_boosted_semileptonic_model(
            ['theta_input_elec_toptag_rebinned.root'],
            narrow_resonances,
            'zp*',
            mcstat,
            eflag = True
        )

    elif type == 'wide_resonances_electron':

        model = build_boosted_semileptonic_model(
            ['theta_input_elec_toptag_rebinned.root'],
            wide_resonances,
            'zp*',
            mcstat,
            eflag = True
        )

    elif type == 'rsg_resonances_electron':

        model = build_boosted_semileptonic_model(
            ['theta_input_elec_toptag_rebinned.root'],
            rsg_resonances,
            'rsg*',
            mcstat,
            eflag = True
        )

    elif type == 'narrow_resonances_lepton':

        model = build_boosted_semileptonic_model(
            ['theta_input_lepton_toptag_rebinned.root'],
            narrow_resonances,
            'zp*',
            mcstat,
            eflag = True,
            muflag = True
        )

    elif type == 'wide_resonances_lepton':

        model = build_boosted_semileptonic_model(
            ['theta_input_lepton_toptag_rebinned.root'],
            wide_resonances,
            'zp*',
            mcstat,
            eflag = True,
            muflag = True
        )

    elif type == 'rsg_resonances_lepton':

        model = build_boosted_semileptonic_model(
            ['theta_input_lepton_toptag_rebinned.root'],
            rsg_resonances,
            'rsg*',
            mcstat,
            eflag = True,
            muflag = True
        )

    else:

        raise exceptions.ValueError('Type %s is undefined' % type)

    for p in model.distribution.get_parameters():
        d = model.distribution.get_distribution(p)
        if d['typ'] == 'gauss' and d['mean'] == 0.0 and d['width'] == 1.0:
            model.distribution.set_distribution_parameters(p, range = [-5.0, 5.0])
        #if 'rate' in p:
        #    if d['typ'] == 'gauss' and d['mean'] == 0.0 and d['width'] == 1.0:
        #        model.distribution.set_distribution_parameters(p, range = [-5.0, 5.0])
        #else:
        #    if d['typ'] == 'gauss' and d['mean'] == 0.0 and d['width'] == 1.0:
        #        model.distribution.set_distribution_parameters(p, range = [-0.0, 0.0])

    return model


# Code introduced by theta_driver

# Building the statistical model
args = {'type': 'narrow_resonances_muon'}

model = build_model(**args)

args = {}

results = bayesian_limits(model, run_theta = True, **args)
exp, obs = results
execfile("utils.py")
limit_table(exp, obs)
