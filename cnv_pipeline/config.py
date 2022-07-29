import os


def _get_gatk_alias():
    """Check for GATK_ALIAS in env, else use `gatk`."""
    gatk_alias = os.environ.get('GATK_ALIAS', 'gatk')
    if '$' in gatk_alias:
        gatk_alias = os.path.expandvars(gatk_alias)
    return gatk_alias


GATK_ALIAS = _get_gatk_alias()
