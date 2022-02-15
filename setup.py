import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="seastatecci_whales",
    version="1.0",
    author="Marcello Passaro",
    author_email="marcello.passaro@tum.de",
    description="WHALES retracker.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitlab.lrz.de/ne62rut/seastatecci_whales.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2.7",
        #"License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        #"Operating System :: OS Independent",
    ],
    python_requires='==2.7*',
    #setup_requires=['wget']
    #entry_points={'console_scripts': ['nctoolbox=nctoolbox.l2_matchup_conversion:__main__'],},
)
