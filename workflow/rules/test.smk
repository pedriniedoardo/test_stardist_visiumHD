# ---------------------------------------------------------------------------- #
rule test:
    '''
    This is a test rule
    '''
    input:
        test = config['test01']
    conda:
        config['env_bin2cell']
    output:
        test = f"{config["out_location"]}{config["test_out"]}"
    log:
        'logs/test/test.log'
    benchmark:
        'benchmarks/test/test.txt'
    resources:
        mem_mb = 500,
        cpus = 1
    threads: 1
    params:
        annotations = config['test02']
    shell:
        '''
        echo "run the test" > {log}
        cat {input.test} > {output.test}
        cat {params.annotations} | awk 'NR<=10' >> {output.test}
        echo "test done" >> {log}
        '''
