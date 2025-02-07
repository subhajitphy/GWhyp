from setuptools import setup, find_packages

setup(
    name="GWhyp",
    version="0.1.1",
    description="Compute Linear GW memory PTA signals from hyperbolic encounters of SMBH system",
    author="Subhajit Dandapat",
    author_email="subhajit.phy97@gmail.com",
    url='https://github.com/subhajitphy/GWhyp',
    python_requires=">=3.8",
    package_dir={"": "./src/"},
    py_modules=["constants","gw_hyp_valid","gw_hyp","gwhyp_v2","hypmik3pn","gw_functions","antenna_pattern","getx","eval_max"],
)
