import setuptools

with open("Readme.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PyFEMP",
    version="0.0.1",
    author="Sascha F. Maassen",
    author_email="sascha.maassen@uni-due.de",
    description="A python FEM solver for educational purposes.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://git.uni-due.de/mechanik/py_femp",
    packages=setuptools.find_packages(),
    project_urls={
        'Wiki' : 'https://git.uni-due.de/mechanik/py_femp/-/wikis/home'
    },
    license='MIT',
    install_requires=[
        "numpy",
        "matplotlib"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)