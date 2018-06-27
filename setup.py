from setuptools import setup

setup(
    version='1.0.0',

    entry_points={
        'console_scripts': [
            'minitv = minitv.__main__:main'
        ]
    }
)
