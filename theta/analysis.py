
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
    model.add_lognormal_uncertainty('wl_rate', math.log(1.5), 'wlight')
    model.add_lognormal_uncertainty('wb_rate', math.log(2.0), 'wb')
    model.add_lognormal_uncertainty('wc_rate', math.log(2.0), 'wc')
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
        #    model.add_lognormal_uncertainty(sys, math.log( _DIFF_ ), proc, obs)
        #  Where _DIFF_ = 2.0 * proc_sys_plus.Integral() / (proc_sys_minus.Integral() + proc_sys_plus.Integral())
        #    after the Chi2 selection
        for obs in ['el_1top0btag_mttbar','el_1top1btag_mttbar','el_1top2btag_mttbar','el_0top0btag_mttbar','el_0top1btag_mttbar','el_0top2btag_mttbar']:
            model.add_lognormal_uncertainty('scale_vjets', math.log(0.50611), 'wb', obs)
            model.add_lognormal_uncertainty('matching_vjets', math.log(0.98630), 'wb', obs)

            model.add_lognormal_uncertainty('scale_vjets', math.log(0.50099), 'wc', obs)
            model.add_lognormal_uncertainty('matching_vjets', math.log(1.01989), 'wc', obs)

        for obs in ['el_1top0btag_mttbar','el_1top1btag_mttbar','el_1top2btag_mttbar','el_0top2btag_mttbar']:
            model.add_lognormal_uncertainty('scale_vjets', math.log(0.49152), 'wlight', obs)
            model.add_lognormal_uncertainty('matching_vjets', math.log(1.02911), 'wlight', obs)

        for p in model.processes:
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
