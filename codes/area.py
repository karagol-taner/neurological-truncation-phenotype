import panel_v2


class InArea:
    pattern = "(%s) and not (%s)" % (panel_v2.DISEASE.pattern, panel_v2.EXCL.pattern)

    def search(self, text):
        return bool(panel_v2.DISEASE.search(text)) and not panel_v2.EXCL.search(text)


IN_AREA = InArea()
