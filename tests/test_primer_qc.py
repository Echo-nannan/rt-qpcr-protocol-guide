from rt_qpcr_guide.primer_qc import PrimerPair, check_primer, check_primer_pair


def test_check_primer_reports_basic_metrics():
    qc = check_primer("ATGCATGCATGCATGCATGC")

    assert qc.length == 20
    assert qc.gc_percent == 50.0
    assert qc.wallace_tm_c == 60
    assert qc.passes_basic_rules is True


def test_check_primer_flags_bad_gc_content():
    qc = check_primer("ATATATATATATATATATAT")

    assert qc.gc_percent == 0
    assert qc.passes_basic_rules is False
    assert "GC content outside 40-60%" in qc.warnings


def test_check_primer_pair_checks_tm_difference_and_product_size():
    pair = PrimerPair(
        forward="ATGCATGCATGCATGCATGC",
        reverse="GCATGCATGCATGCATGCAT",
        product_size_bp=120,
    )

    qc = check_primer_pair(pair)

    assert qc.product_size_ok is True
    assert qc.tm_difference_ok is True
    assert qc.passes is True
