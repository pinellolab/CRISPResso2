import logging
from CRISPResso2.CRISPRessoShared import (
    GuardrailMessageHandler,
    TotalReadsGuardrail,
    OverallReadsAlignedGuardrail,
    DisproportionateReadsAlignedGuardrail,
    LowRatioOfModsInWindowToOutGuardrail,
    HighRateOfModificationAtEndsGuardrail,
    HighRateOfSubstitutionsOutsideWindowGuardrail,
    HighRateOfSubstitutionsGuardrail,
    ShortSequenceGuardrail,
    LongAmpliconShortReadsGuardrail,
)


def make_handler():
    logger = logging.getLogger("test_guardrails")
    return GuardrailMessageHandler(logger)


# --- GuardrailMessageHandler ---

class TestGuardrailMessageHandler:
    def test_initial_state(self):
        handler = make_handler()
        assert handler.get_failed_html_messages() == []
        assert handler.get_passed_html_messages() == []
        assert handler.get_messages() == {}

    def test_report_warning_creates_danger_alert(self):
        handler = make_handler()
        handler.report_warning("Something went wrong.")
        msgs = handler.get_failed_html_messages()
        assert len(msgs) == 1
        assert "alert-danger" in msgs[0]
        assert "Something went wrong." in msgs[0]

    def test_report_pass_creates_success_alert(self):
        handler = make_handler()
        handler.report_pass("All good.")
        msgs = handler.get_passed_html_messages()
        assert len(msgs) == 1
        assert "alert-success" in msgs[0]
        assert "All good." in msgs[0]

    def test_display_warning_stores_message(self):
        handler = make_handler()
        handler.display_warning("TestGuardrail", "a warning")
        assert handler.get_messages() == {"TestGuardrail": "a warning"}

    def test_warning_html_has_space_after_strong(self):
        handler = make_handler()
        handler.report_warning("Test message.")
        html = handler.get_failed_html_messages()[0]
        assert "</strong> Test message." in html


# --- TotalReadsGuardrail ---

class TestTotalReadsGuardrail:
    def test_fail_below_minimum(self):
        handler = make_handler()
        g = TotalReadsGuardrail(handler, 10000)
        g.safety(500)
        assert len(handler.get_failed_html_messages()) == 1
        assert len(handler.get_passed_html_messages()) == 0
        assert "500" in handler.get_failed_html_messages()[0]

    def test_pass_above_minimum(self):
        handler = make_handler()
        g = TotalReadsGuardrail(handler, 10000)
        g.safety(20000)
        assert len(handler.get_failed_html_messages()) == 0
        assert len(handler.get_passed_html_messages()) == 1

    def test_fail_at_exact_minimum(self):
        handler = make_handler()
        g = TotalReadsGuardrail(handler, 10000)
        g.safety(10000)
        # 10000 is not < 10000, so it should pass
        assert len(handler.get_failed_html_messages()) == 0
        assert len(handler.get_passed_html_messages()) == 1


# --- OverallReadsAlignedGuardrail ---

class TestOverallReadsAlignedGuardrail:
    def test_fail_low_alignment(self):
        handler = make_handler()
        g = OverallReadsAlignedGuardrail(handler, 0.9)
        g.safety(1000, 500)  # 50% aligned
        assert len(handler.get_failed_html_messages()) == 1
        assert len(handler.get_passed_html_messages()) == 0

    def test_pass_high_alignment(self):
        handler = make_handler()
        g = OverallReadsAlignedGuardrail(handler, 0.9)
        g.safety(1000, 950)  # 95% aligned
        assert len(handler.get_failed_html_messages()) == 0
        assert len(handler.get_passed_html_messages()) == 1

    def test_zero_total_reads(self):
        handler = make_handler()
        g = OverallReadsAlignedGuardrail(handler, 0.9)
        g.safety(0, 0)
        # Should return early, no messages
        assert len(handler.get_failed_html_messages()) == 0
        assert len(handler.get_passed_html_messages()) == 0

    def test_fail_at_exact_cutoff(self):
        handler = make_handler()
        g = OverallReadsAlignedGuardrail(handler, 0.9)
        g.safety(1000, 900)  # exactly 90%
        # <= cutoff means fail
        assert len(handler.get_failed_html_messages()) == 1


# --- DisproportionateReadsAlignedGuardrail ---

class TestDisproportionateReadsAlignedGuardrail:
    def test_single_amplicon_skipped(self):
        handler = make_handler()
        g = DisproportionateReadsAlignedGuardrail(handler, 0.3)
        g.safety(1000, {"amp1": 1000})
        # Single amplicon should be skipped entirely
        assert len(handler.get_failed_html_messages()) == 0
        assert len(handler.get_passed_html_messages()) == 0

    def test_even_distribution_triggers_warning(self):
        # NOTE: This documents a known bug in the upper-bound check.
        # The condition `aligned >= expected * (1 - cutoff)` with cutoff=0.9
        # becomes `aligned >= expected * 0.1`, which nearly always triggers.
        # An even 50/50 split incorrectly fails.
        handler = make_handler()
        g = DisproportionateReadsAlignedGuardrail(handler, 0.9)
        g.safety(1000, {"amp1": 500, "amp2": 500})
        assert len(handler.get_failed_html_messages()) == 2
        assert len(handler.get_passed_html_messages()) == 0

    def test_fail_heavily_skewed(self):
        handler = make_handler()
        g = DisproportionateReadsAlignedGuardrail(handler, 0.3)
        g.safety(1000, {"amp1": 950, "amp2": 50})
        # amp2 has 50 reads, expected is 500, cutoff 0.3 => lower bound 150
        # 50 <= 150, so amp2 fails
        assert len(handler.get_failed_html_messages()) >= 1

    def test_pass_balanced_low_cutoff(self):
        # With a very low cutoff (0.01), the upper bound check becomes
        # aligned >= expected * 0.99, which a perfectly even split still exceeds.
        # This further documents the bug — even cutoff=0.01 triggers on even splits.
        handler = make_handler()
        g = DisproportionateReadsAlignedGuardrail(handler, 0.01)
        g.safety(1000, {"amp1": 500, "amp2": 500})
        assert len(handler.get_failed_html_messages()) == 2


# --- LowRatioOfModsInWindowToOutGuardrail ---

class TestLowRatioOfModsInWindowToOutGuardrail:
    def test_fail_low_ratio(self):
        handler = make_handler()
        g = LowRatioOfModsInWindowToOutGuardrail(handler, 0.01)
        g.safety(0, 1000)  # 0% in window
        assert len(handler.get_failed_html_messages()) == 1
        assert len(handler.get_passed_html_messages()) == 0

    def test_pass_high_ratio(self):
        handler = make_handler()
        g = LowRatioOfModsInWindowToOutGuardrail(handler, 0.01)
        g.safety(900, 100)  # 90% in window
        assert len(handler.get_failed_html_messages()) == 0
        assert len(handler.get_passed_html_messages()) == 1

    def test_zero_total_mods(self):
        handler = make_handler()
        g = LowRatioOfModsInWindowToOutGuardrail(handler, 0.01)
        g.safety(0, 0)
        assert len(handler.get_failed_html_messages()) == 0
        assert len(handler.get_passed_html_messages()) == 0


# --- HighRateOfModificationAtEndsGuardrail ---

class TestHighRateOfModificationAtEndsGuardrail:
    def test_fail_high_rate(self):
        handler = make_handler()
        g = HighRateOfModificationAtEndsGuardrail(handler, 0.01)
        g.safety(1000, 100)  # 10% irregular
        assert len(handler.get_failed_html_messages()) == 1
        assert len(handler.get_passed_html_messages()) == 0

    def test_pass_low_rate(self):
        handler = make_handler()
        g = HighRateOfModificationAtEndsGuardrail(handler, 0.01)
        g.safety(1000, 5)  # 0.5% irregular
        assert len(handler.get_failed_html_messages()) == 0
        assert len(handler.get_passed_html_messages()) == 1

    def test_zero_total_reads(self):
        handler = make_handler()
        g = HighRateOfModificationAtEndsGuardrail(handler, 0.01)
        g.safety(0, 0)
        assert len(handler.get_failed_html_messages()) == 0
        assert len(handler.get_passed_html_messages()) == 0


# --- HighRateOfSubstitutionsOutsideWindowGuardrail ---

class TestHighRateOfSubstitutionsOutsideWindowGuardrail:
    def test_fail_high_rate(self):
        handler = make_handler()
        g = HighRateOfSubstitutionsOutsideWindowGuardrail(handler, 0.002)
        g.safety(1000, 500)  # 50% outside
        assert len(handler.get_failed_html_messages()) == 1
        assert len(handler.get_passed_html_messages()) == 0

    def test_pass_low_rate(self):
        handler = make_handler()
        g = HighRateOfSubstitutionsOutsideWindowGuardrail(handler, 0.002)
        g.safety(1000, 1)  # 0.1% outside
        assert len(handler.get_failed_html_messages()) == 0
        assert len(handler.get_passed_html_messages()) == 1

    def test_zero_global_subs(self):
        handler = make_handler()
        g = HighRateOfSubstitutionsOutsideWindowGuardrail(handler, 0.002)
        g.safety(0, 0)
        assert len(handler.get_failed_html_messages()) == 0
        assert len(handler.get_passed_html_messages()) == 0


# --- HighRateOfSubstitutionsGuardrail ---

class TestHighRateOfSubstitutionsGuardrail:
    def test_fail_high_sub_rate(self):
        handler = make_handler()
        g = HighRateOfSubstitutionsGuardrail(handler, 0.3)
        g.safety(500, 500, 400)  # 400/1000 = 40% subs
        assert len(handler.get_failed_html_messages()) == 1
        assert len(handler.get_passed_html_messages()) == 0

    def test_pass_low_sub_rate(self):
        handler = make_handler()
        g = HighRateOfSubstitutionsGuardrail(handler, 0.3)
        g.safety(500, 500, 100)  # 100/1000 = 10% subs
        assert len(handler.get_failed_html_messages()) == 0
        assert len(handler.get_passed_html_messages()) == 1

    def test_zero_total_mods(self):
        handler = make_handler()
        g = HighRateOfSubstitutionsGuardrail(handler, 0.3)
        g.safety(0, 0, 0)
        assert len(handler.get_failed_html_messages()) == 0
        assert len(handler.get_passed_html_messages()) == 0


# --- ShortSequenceGuardrail ---

class TestShortSequenceGuardrail:
    def test_fail_short_amplicon(self):
        handler = make_handler()
        g = ShortSequenceGuardrail(handler, 50, "amplicon")
        g.safety({"amp1": 30})
        assert len(handler.get_failed_html_messages()) == 1
        assert "amp1" in handler.get_failed_html_messages()[0]

    def test_pass_long_amplicon(self):
        handler = make_handler()
        g = ShortSequenceGuardrail(handler, 50, "amplicon")
        g.safety({"amp1": 200})
        assert len(handler.get_failed_html_messages()) == 0
        assert len(handler.get_passed_html_messages()) == 1

    def test_fail_short_guide(self):
        handler = make_handler()
        g = ShortSequenceGuardrail(handler, 19, "guide")
        g.safety({"ATCG": 4})
        assert len(handler.get_failed_html_messages()) == 1
        assert len(handler.get_passed_html_messages()) == 0

    def test_pass_long_guide(self):
        handler = make_handler()
        g = ShortSequenceGuardrail(handler, 19, "guide")
        g.safety({"ATCGATCGATCGATCGATCG": 20})
        assert len(handler.get_failed_html_messages()) == 0
        assert len(handler.get_passed_html_messages()) == 1

    def test_multiple_sequences_mixed(self):
        handler = make_handler()
        g = ShortSequenceGuardrail(handler, 50, "amplicon")
        g.safety({"amp1": 200, "amp2": 30, "amp3": 100})
        # One failure (amp2), but one pass message overall only if all pass
        assert len(handler.get_failed_html_messages()) == 1
        assert len(handler.get_passed_html_messages()) == 0

    def test_message_does_not_start_with_space(self):
        handler = make_handler()
        g = ShortSequenceGuardrail(handler, 50, "amplicon")
        g.safety({"amp1": 30})
        html = handler.get_failed_html_messages()[0]
        # Extract the message after the strong tag
        assert "</strong> Amplicon" in html


# --- LongAmpliconShortReadsGuardrail ---

class TestLongAmpliconShortReadsGuardrail:
    def test_fail_long_amplicon(self):
        handler = make_handler()
        g = LongAmpliconShortReadsGuardrail(handler, 1.5)
        g.safety({"amp1": 300}, 100)  # 300 > 100*1.5
        assert len(handler.get_failed_html_messages()) == 1
        assert len(handler.get_passed_html_messages()) == 0

    def test_pass_matching_lengths(self):
        handler = make_handler()
        g = LongAmpliconShortReadsGuardrail(handler, 1.5)
        g.safety({"amp1": 140}, 100)  # 140 < 100*1.5
        assert len(handler.get_failed_html_messages()) == 0
        assert len(handler.get_passed_html_messages()) == 1

    def test_multiple_amplicons_one_fails(self):
        handler = make_handler()
        g = LongAmpliconShortReadsGuardrail(handler, 1.5)
        g.safety({"amp1": 100, "amp2": 500}, 100)
        assert len(handler.get_failed_html_messages()) == 1
        assert "amp2" in handler.get_failed_html_messages()[0]
        assert len(handler.get_passed_html_messages()) == 0
